import Pkg
Pkg.activate(".")
using RCall
begin
    R"""
    library(conflicted)
    library(tidyverse)
    library(here)
    library(embarcadero)
    library(viridis)
    library(furrr)
    plan(multiprocess)

    # Resolve conflicts
    conflict_prefer("filter", "dplyr")
    conflict_prefer("intersect", "dplyr")
    conflict_prefer("select", "dplyr")

    # Custom functions
    source(here("src", "lib", "R", "bart.R"))
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

    message("Loading & preparing data")

    # Load data
    spa_full <- read_tsv(here("data", "proc", "distributions_spa_full.csv"))
    env_full <- read_tsv(here("data", "proc", "distributions_env_full.csv"))
    spe      <- read_tsv(here("data", "proc", "distributions_spe_full.csv"))

    # Load QC data (optional)
    spa_qc <- read_tsv(here("data", "proc", "distributions_spa_qc.csv"))

    # Prepare data
    # subset_qc <- TRUE # subset to QC data (optional)
    source(here("src", "02_training_data-preparation.R"))

    # Remove site with NAs for landcover variables
    if (length(inds_withNAs) > 0) {
        richness_values <- richness_values[-inds_withNAs]
        lcbd_values <- lcbd_values[-inds_withNAs]
    }

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


# Map richness difference
richness_diff = similar(richness_raw)
richness_diff.grid = abs.(richness_bart.grid .- richness_raw.grid)
diff_plot = plotSDM(richness_diff, c = :inferno, clim = (-Inf, Inf),
                    title = "Predicted richness - BART vs raw",
                    colorbar_title = "Difference in predicted richness (absolute)",
                    )
histogram(filter(!isnan, richness_diff.grid), bins = 20)

## Predictions for full range
begin
    R"""
    # Remove sites with NA values
    inds_na <- map(env_full, ~ which(is.na(.x)))
    (inds_na <- sort(unique(unlist(inds_na))))
    vars_nona <- vars_full[-inds_na,]

    # Make predictions
    predictions_full_scale <- predict(classif_model, vars_nona)$predictions

    # Add sites with NAs
    predictions_full <- matrix(NA, nrow = nrow(vars_full), ncol = 1)
    colnames(predictions_full) <- colnames(predictions_full_scale)
    predictions_full[-inds_na,] <- predictions_full_scale

    """
end

@rget predictions_full

# Arrange as layer
replace!(predictions_full, missing => NaN)
predictions_full = reshape(predictions_full, size(richness_raw.grid)) |> Array{Float64}
richness_rf_full = similar(richness_raw)
richness_rf_full.grid = predictions_full

# Visualize result
richness_plot_full = plotSDM(richness_rf_full, c = :viridis,
                             title = "Richness RF predictions- All sites",
                             colorbar_title = "Predicted number of species",
                             )

# Get comparison
@load joinpath("data/", "jld2", "rf-distributions.jld2") distributions
Ysdm = calculate_Y(distributions)
richness_sdm = calculate_richness(Ysdm, distributions[1])
plotSDM(richness_sdm, c = :viridis)

# Map richness difference
richness_diff_full = similar(richness_rf_full)
richness_diff_full.grid = abs.(richness_rf_full.grid .- richness_sdm.grid)
diff_plot_full = plotSDM(richness_diff_full, c = :inferno, clim = (-Inf, Inf),
                         title = "Predicted richness - RF vs SDM",
                         colorbar_title = "Difference in predicted richness (absolute)",
                         )
histogram(filter(!isnan, richness_diff.grid), bins = 20)

## Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(richness_plot, dpi = 150),      joinpath("fig", "raw", "x_raw_richness-rf.png"))
    savefig(plot(richness_plot_full, dpi = 150), joinpath("fig", "rf",  "x_rf_richness-rf.png"))

    savefig(plot(diff_plot, dpi = 150),      joinpath("fig", "raw", "x_raw_richness-diff.png"))
    savefig(plot(diff_plot_full, dpi = 150), joinpath("fig", "rf",  "x_rf_richness-diff.png"))
end
