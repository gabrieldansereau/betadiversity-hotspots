import Pkg
Pkg.activate(".")
using RCall
begin
    R"""
    source(file.path("src", "required.R"))
    """
end
using Distributed
@time include(joinpath("..", "required.jl"))

## Conditional arguments
# save_figures = true
# subset_qc = true

## Load distribution data for all species
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions

# Subset to QC data (optional)
coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
if (@isdefined subset_qc)
    distributions = [d[coords_qc] for d in distributions]
    @rput subset_qc
end

## Richness
Y = calculate_Y(distributions)
richness_raw = calculate_richness(Y, distributions[1])
lcbd_raw = calculate_lcbd(Y, distributions[1])

# Visualize
plotSDM(richness_raw, c = :viridis)
plotSDM(lcbd_raw, c = :viridis)
histogram2d(richness_raw, lcbd_raw, c = :viridis, bin = 40, ylim = (0.0, 1.0))

## Extract values for model
inds_obs = _indsobs(Y)
richness_values = Int64.(richness_raw.grid[inds_obs])
lcbd_values = lcbd_raw.grid[inds_obs]

## Train BART
@rput richness_values lcbd_values inds_obs
begin
    R"""
    ## 1. Load data ####

    # Load data
    source(here("src", "02a_training_data-preparation.R"))

    # Remove site with NAs for landcover variables
    if (length(inds_withNAs) > 0) {
        lcbd_values <- lcbd_values[-inds_withNAs]
    }
    
    # Get richness values
    richness_values <- rowSums(spe)
    
    # Get LCBD values
    # Apply Hellinger transformation
    spe_transf <- vegan::decostand(spe_full[inds_obs,-1], "hel")
    # S -> squared deviations from column mean
    S <- map_df(spe_transf, ~ (.x - mean(.x))^2)
    # SStotal -> total sum of squares
    SStotal <- sum(S)
    # SSi -> sum of squares for site i
    SSi <- rowSums(S)
    # LCBD -> local contribution to beta diversity (site i, relative)
    lcbd_values_R <- SSi/SStotal
    lcbd_values_R <- lcbd_values_R/max(lcbd_values_R)
    
    if (length(inds_withNAs) > 0) {
        lcbd_values_R <- lcbd_values_R[-inds_withNAs]
    }

    length(lcbd_values_R) == length(lcbd_values)
    lcbd_values_R == lcbd_values
    all(near(lcbd_values_R, lcbd_values))
    sum(round(lcbd_values_R, 5) == round(lcbd_values, 5))
    

    # Select variables
    xnames <- select(vars, -c(lat, lon)) %>% names()

    ## 2. Create layers ####

    message("Creating layers")

    # Create raster layers
    vars_layers <- map(
        vars_full[,xnames], 
        ~ df_to_layer(.x, lons = vars_full$lon, lats = vars_full$lat)
    )
    wc_layer <- vars_layers$wc1

    # Stack variables layers
    vars_stack <- stack(vars_layers, names = xnames)


    ## 3. Model training ####

    # Group richness & lcbd
    values_df <- tibble(
        richness = richness_values,
        lcbd = lcbd_values
    )
    
    # Function to run BART in parallel
    bart_parallel <- function(...) {
        # BART model
        model <- bart(...)
        # Touch state so that saving will work
        invisible(model$fit$state)
        return(model)
    }

    # Train model
    set.seed(42)
    system.time(
        models <- future_map(
            values_df,
            ~ bart_parallel(
                y.train = .x,
                x.train = vars[,xnames],
                keeptrees = TRUE
            )
        )
    ) # 90 sec.
    varimp(models[[1]], plot = TRUE)
    varimps <-  map_dfr(models, varimp, .id = "value") %>% 
        rename(vars = names)
    varimps %>% 
        pivot_wider(names_from = "value", values_from = "varimps") %>% 
        mutate(diff = richness - lcbd)
    ggplot(varimps, aes(vars)) + 
        geom_bar(aes(weight = varimps, fill = value)) #, position = "dodge2")

    ggplot(varimps, aes(vars, varimps)) + 
        geom_boxplot(aes(colour = value))

    ## 4. Predictions ####
    system.time(
        predictions <- future_map(
            models,
            ~ predict2.bart(
                .x,
                vars_stack,
                quantiles = c(0.025, 0.0975),
                splitby = 20
            )
        )
    ) # 2 min.

    # Plot richness
    plot(
        predictions$richness[[1]], 
        main = 'Probability predictions - Posterior mean', 
        col = viridis(255),
        # zlim = c(0, 1),
        legend.args=list(text='Probability', side=2, line=1.3),
        box = FALSE, axes = FALSE
    )
    # Negative richness??
    plot(
        predictions$richness[[1]] < 1, 
        main = 'Probability predictions - Posterior mean', 
        col = viridis(255),
        # zlim = c(0, 1),
        legend.args=list(text='Probability', side=2, line=1.3),
        box = FALSE, axes = FALSE
    )
    plot(
        predictions$lcbd[[1]], 
        main = 'Probability predictions - Posterior mean', 
        col = viridis(255),
        zlim = c(0, 1),
        legend.args=list(text='Probability', side=2, line=1.3),
        box = FALSE, axes = FALSE
    )

    # Predictions
    pred_df <- predictions %>% 
        map(~ .x$layer.1) %>% 
        stack() %>% 
        as.data.frame(xy = TRUE) %>% 
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    pred_df
    # Lower quantiles
    lower_df <- predictions %>% 
        map(~ .x$layer.2) %>% 
        stack() %>% 
        as.data.frame(xy = TRUE) %>% 
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    lower_df
    # Upper quantiles
    upper_df <- predictions %>% 
        map(~ .x$layer.3) %>%
        stack() %>% 
        as.data.frame(xy = TRUE) %>%  
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    upper_df
    """
end
@rget pred_df lower_df upper_df

# Fix missing values
predictions = replace(Array(pred_df), missing => NaN)
lower = replace(Array(pred_df), missing => NaN)
upper = replace(Array(pred_df), missing => NaN)

## Plot predicted richness
richness_bart = similar(richness_raw)
richness_bart.grid[:] = predictions[:, 1]
richness_plot = plotSDM(richness_bart, c = :viridis,
                        title = "Richness BART predictions",
                        colorbar_title = "Predicted number of species",
                        )

## Plot predicted LCBD
lcbd_bart = similar(lcbd_raw)
lcbd_bart.grid[:] = predictions[:, 2]
lcbd_plot = plotSDM(lcbd_bart, c = :viridis,
                    title = "LCBD BART predictions",
                    colorbar_title = "LCBD scores",
                    clim = (0,1),
                    )

## Get full-scale comparison
@load joinpath("data/", "jld2", "bart-distributions.jld2") distributions
Ysdm = calculate_Y(distributions)
richness_sdm = calculate_richness(Ysdm, distributions[1])
lcbd_sdm = calculate_lcbd(Ysdm, distributions[1])

# Map richness difference
richness_diff = similar(richness_bart)
richness_diff.grid = abs.(richness_bart.grid .- richness_sdm.grid)
richness_diff_plot = plotSDM(richness_diff, c = :inferno, clim = (-Inf, Inf),
                         title = "Predicted richness - BART vs SDM",
                         colorbar_title = "Difference in predicted richness (absolute)",
                         )
histogram(filter(!isnan, richness_diff.grid), bins = 20)

# Map LCBD difference
lcbd_diff = similar(lcbd_bart)
lcbd_diff.grid = abs.(lcbd_bart.grid .- lcbd_sdm.grid)
lcbd_diff_plot = plotSDM(lcbd_diff, c = :inferno, clim = (-Inf, Inf),
                         title = "Predicted LCBD - BART vs SDM",
                         colorbar_title = "Difference in predicted LCBD (absolute)",
                         )
histogram(filter(!isnan, lcbd_diff.grid), bins = 20)

## Export figures
save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(richness_plot, dpi = 150), joinpath("fig", "bart", "x_bart_richness-bart.png"))
    savefig(plot(lcbd_plot, dpi = 150),     joinpath("fig", "bart", "x_bart_lcbd-bart.png"))

    savefig(plot(richness_diff_plot, dpi = 150), joinpath("fig", "bart", "x_bart_richness-diff.png"))
    savefig(plot(lcbd_diff_plot, dpi = 150),     joinpath("fig", "bart", "x_bart_lcbd-diff.png"))
end
