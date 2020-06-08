#### 03b - Random Forests predictions ####
import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

## Conditional arguments
# save_data = true

## Perform RandomForests
using RCall
begin
    R"""
    ## 0. Load packages ####
    source(file.path("src", "required.R"))

    # Conditional evaluations
    # subset_qc <- TRUE # subset to QC data (optional)
    # create_models <- TRUE # train models
    # save_models <- TRUE # save & overwrite models

    ## 1. Load data ####

    # Prepare data
    source(here("src", "02a_training_data-preparation.R"))

    # Remove sites with NA values
    inds_na <- map(env_full, ~ which(is.na(.x)))
    (inds_na <- sort(unique(unlist(inds_na))))
    vars_nona <- vars_full[-inds_na,]

    # Load model
    load(here("data", "proc", "rf_models.RData"))

    # Make prediction
    system.time(
        rf_pred <- pblapply(ranger_models, function(x) predict(x, vars_nona))
    )

    # Extract predictions
    predictions <- map(rf_pred, "predictions") %>% 
        map_df(~ as.numeric(levels(.x))[.x])

    # Add sites with NAs
    predictions_full <- matrix(
        NA, 
        nrow = nrow(vars_full), 
        ncol = ncol(predictions)
    )
    colnames(predictions_full) <- colnames(predictions)
    predictions_full[-inds_na,] <- as.matrix(predictions)
    """
end
@rget predictions predictions_full inds_na

## Create Y matrices

# Get matrix Y
Y = replace(predictions_full, missing => NaN)
# Set values to NaN if no species present
inds_zeros = findall(map(x -> all(iszero.(x)), eachrow(Y)))
Y[inds_zeros,:] .= NaN

## Create distributions

# Load raw distributions (for grid size)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
raw_distributions = distributions
# Get layer dimensions & limits
dims = size(raw_distributions[1].grid)
lims = (left = raw_distributions[1].left, right = raw_distributions[1].right,
        bottom = raw_distributions[1].bottom, top = raw_distributions[1].top)

# Create RF distribution layers
Ydistrib = replace(Y, 0.0 => NaN)
Ygrids = [Ydistrib[:, col] for col in 1:size(Ydistrib,2)]
Ygrids = reshape.(Ygrids, dims...)
distributions = SimpleSDMResponse.(Ygrids, lims...)

## Export results
# save_data = true
if @isdefined save_data && save_data == true
    @save joinpath("data", "jld2", "rf-distributions.jld2") distributions
    _zip_jld2(joinpath("data", "jld2", "rf-distributions.zip"),
              joinpath("data", "jld2", "rf-distributions.jld2"))
    touch(joinpath("data", "jld2", "rf-distributions.jld2"))
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])
