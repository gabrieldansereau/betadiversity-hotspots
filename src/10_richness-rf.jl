import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# save_figures = true

## Load distribution data for all species
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions spenames speindex

## Load matrix Y
@load joinpath("data", "jld2", "raw-Y-matrices.jld2") Y Yobs Ytransf inds_obs inds_notobs

## Richness
richness_raw = calculate_richness(Y, inds_notobs, distributions)

# Visualize
plotSDM(richness_raw, c = :viridis)

## Extract values for model
richness_values = Int64.(richness_raw.grid[inds_obs])

## Train Random Forest
@rput richness_values inds_obs
begin
    R"""
    library(ranger)
    spa <- read.csv("data/proc/distributions_spa.csv", header=TRUE, sep="\t")
    env <- read.csv("data/proc/distributions_env.csv", header=TRUE, sep="\t")

    # Remove site with NAs for landcover variables
    (inds_withNAs <- unique(unlist(sapply(env, function(x) which(is.na(x))))))
    if (length(inds_withNAs) > 0) {
        richness_values <- richness_values[-inds_withNAs]
        spa <- spa[-inds_withNAs,]
        env <- env[-inds_withNAs,]
    }

    # Combine environmental variables
    vars <- cbind(env, spa)

    # Separate into training/testing datasets
    set.seed(42)
    inds_train <- sample(nrow(vars), 0.7*nrow(vars), replace = FALSE)

    richness_train <- richness_values[inds_train]
    vars_train <- vars[inds_train,]

    richness_test <- richness_values[-inds_train]
    vars_test <- vars[-inds_train,]

    # Train model
    system.time(regress_model <- ranger(richness_train ~ ., data = vars_train, importance = "impurity", seed = 42))
    system.time(classif_model <- ranger(as.factor(richness_train) ~ ., data = vars_train, importance = "impurity", seed = 42))

    regress_pred <- predict(regress_model, vars_test)$predictions
    sum(round(regress_pred) == richness_test)/length(richness_test)

    classif_pred <- predict(classif_model, vars_test)$predictions
    sum(classif_pred == richness_test)/length(richness_test)

    """
end

## Extract Random Forest predictions for observed sites
begin
    R"""
    predictions <- predict(classif_model, vars)$predictions
    predictions <- as.numeric(levels(predictions))[predictions]
    """
end
@rget predictions inds_withNAs

## Plot predicted richness
richness_rf = similar(richness_raw)
richness_rf.grid[inds_obs[Not(inds_withNAs)]] = predictions

richness_plot = plotSDM(richness_rf, c = :viridis,
                        title = "Richness RF predictions - Observed sites",
                        colorbar_title = "Predicted number of species",
                        dpi = 150)

# Map richness difference
richness_diff = similar(richness_raw)
richness_diff.grid = abs.(richness_rf.grid .- richness_raw.grid)
diff_plot = plotSDM(richness_diff, c = :inferno, clim = (-Inf, Inf),
                    title = "Predicted richness - RF vs raw",
                    colorbar_title = "Difference in predicted richness (absolute)",
                    dpi = 150)
histogram(filter(!isnan, richness_diff.grid), bins = 20)

## Predictions for full range
begin
    R"""
    ## Predict distributions for full range
    env_full <- read.csv("data/proc/distributions_env_full.csv", header = TRUE, sep = "\t")
    spa_full <- read.csv("data/proc/distributions_spa_full.csv", header = TRUE, sep = "\t")
    vars_full <- cbind(env_full, spa_full)
    head(vars_full)

    # Remove sites with NA values
    inds_na <- sapply(env_full, function(x) which(is.na(x)))
    (inds_na <- sort(unique(unlist(inds_na))))
    vars_nona <- vars_full[-inds_na,]

    # Make predictions
    predictions <- predict(classif_model, vars_nona)$predictions
    predictions <- as.numeric(levels(predictions))[predictions]

    # Add sites with NAs
    predictions_full <- matrix(NA, nrow = nrow(vars_full), ncol = 1)
    colnames(predictions_full) <- colnames(predictions)
    predictions_full[-inds_na,] <- predictions

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
                             dpi = 150)

# Get comparison
@load joinpath("data/", "jld2", "rf-distributions.jld2") distributions
sdm = calculate_Ymatrix(distributions)
richness_sdm = calculate_richness(sdm.Y, sdm.inds_notobs, distributions)
plotSDM(richness_sdm, c = :viridis)

# Map richness difference
richness_diff_full = similar(richness_rf_full)
richness_diff_full.grid = abs.(richness_rf_full.grid .- richness_sdm.grid)
diff_plot_full = plotSDM(richness_diff_full, c = :inferno, clim = (-Inf, Inf),
                         title = "Predicted richness - RF vs SDM",
                         colorbar_title = "Difference in predicted richness (absolute)",
                         dpi = 150)
histogram(filter(!isnan, richness_diff.grid), bins = 20)

## Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(richness_plot, joinpath("fig", "raw", "10_raw_richness-rf.png"))
    savefig(richness_plot_full, joinpath("fig", "rf", "10_rf_richness-rf.png"))

    savefig(diff_plot, joinpath("fig", "raw", "10_raw_richness-diff.png"))
    savefig(diff_plot_full, joinpath("fig", "rf", "10_rf_richness-diff.png"))
end
