#### 03b - Random Forests predictions ####
import Pkg; Pkg.activate(".")
using RCall
R"source(file.path('src', 'required.R'))" # bug with `velox` if not called here
include("required.jl")

## Conditional arguments
# save_data = true

## Perform RandomForests
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
    load(here("data", "rdata", "rf_models.RData"))

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
Y = replace(predictions_full, missing => nothing) |> Array{Union{Nothing, Float32}}
# Set values to nothing if no species present
inds_zeros = _indsnotobs(Y)
Y[inds_zeros, :] .= nothing

## Create distributions

# Load raw distributions (for grid size)
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
spenames = filter(:type => ==("species"), glossary).full_name
distributions = [geotiff(SimpleSDMPredictor, joinpath("data", "proc", "distributions_raw.tif"), i) for i in eachindex(spenames)]
raw_distributions = distributions
# Get layer dimensions & limits
dims = size(raw_distributions[1].grid)
lims = (left = raw_distributions[1].left, right = raw_distributions[1].right,
        bottom = raw_distributions[1].bottom, top = raw_distributions[1].top)

# Create RF distribution layers
Ydistrib = replace(Y, 0.0 => nothing)
Ygrids = [Ydistrib[:, col] for col in 1:size(Ydistrib,2)]
Ygrids = reshape.(Ygrids, dims...) .|> Array
distributions = SimpleSDMResponse.(Ygrids, lims...)

## Export results
# save_data = true
if (@isdefined save_data) && save_data == true
    @save joinpath("data", "jld2", "rf-distributions.jld2") distributions
    _zip_jld2(joinpath("data", "jld2", "rf-distributions.zip"),
              joinpath("data", "jld2", "rf-distributions.jld2"))
    touch(joinpath("data", "jld2", "rf-distributions.jld2"))
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])
