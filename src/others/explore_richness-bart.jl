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
subset_qc = true

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

# Visualize
plotSDM(richness_raw, c = :viridis)

## Extract values for model
inds_obs = _indsobs(Y)
richness_values = Int64.(richness_raw.grid[inds_obs])

## Train BART
@rput richness_values inds_obs
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

    # Train model
    set.seed(42)
    system.time(
        model <- bart(
            y.train = richness_values,
            x.train = vars[,xnames],
            keeptrees = TRUE
        )
    ) # 90 sec.
    invisible(model$fit$state)
    summary(model)
    varimp(model, plot = TRUE)


    ## 4. Predictions ####
    system.time(
        predictions <- predict2.bart(
            model,
            vars_stack,
            quantiles = c(0.025, 0.0975),
            splitby = 20
        )
    ) # 2 min.

    # Plot results
    plot(
        predictions[[1]], 
        main = 'Probability predictions - Posterior mean', 
        col = viridis(255),
        # zlim = c(0, 1),
        legend.args=list(text='Probability', side=2, line=1.3),
        box = FALSE, axes = FALSE
    )
    
    plot(
        predictions[[1]] < 1, 
        main = 'Probability predictions - Posterior mean', 
        col = viridis(255),
        # zlim = c(0, 1),
        legend.args=list(text='Probability', side=2, line=1.3),
        box = FALSE, axes = FALSE
    )
    """
end

## Extract Random Forest predictions for observed sites
begin
    R"""
    predictions <- predict(classif_model, vars)$predictions
    """
end
@rget predictions inds_withNAs

## Plot predicted richness
richness_rf = similar(richness_raw)
richness_rf.grid[inds_obs[Not(inds_withNAs)]] = predictions

richness_plot = plotSDM(richness_rf, c = :viridis,
                        title = "Richness RF predictions - Observed sites",
                        colorbar_title = "Predicted number of species",
                        )

# Map richness difference
richness_diff = similar(richness_raw)
richness_diff.grid = abs.(richness_rf.grid .- richness_raw.grid)
diff_plot = plotSDM(richness_diff, c = :inferno, clim = (-Inf, Inf),
                    title = "Predicted richness - RF vs raw",
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
