import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# save_figures = true

## Load data
#=
spe = CSV.read(joinpath("data", "proc", "distributions_spe.csv"), header=true, delim="\t")
spa = CSV.read(joinpath("data", "proc", "distributions_spa.csv"), header=true, delim="\t")
env = CSV.read(joinpath("data", "proc", "distributions_env.csv"), header=true, delim="\t")
=#

## Perform RandomForests
using RCall
# @rput spe spa env
begin
    R"""
    library(ranger)
    library(pbapply)

    ## Predict distributions for full range
    env_full <- read.csv("data/proc/distributions_env_full.csv", header = TRUE, sep = "\t")
    spa_full <- read.csv("data/proc/distributions_spa_full.csv", header = TRUE, sep = "\t")
    vars_full <- cbind(env_full, spa_full)
    head(vars_full)

    # Remove sites with NA values
    inds_na <- sapply(env_full, function(x) which(is.na(x)))
    (inds_na <- sort(unique(unlist(inds_na))))
    vars_nona <- vars_full[-inds_na,]

    # Load model
    load("data/proc/rf_models.RData")

    # Make predictions
    system.time(rf_pred <- pblapply(ranger_models, function(x) predict(x, vars_nona)))

    # Extract predictions
    predictions <- sapply(rf_pred, function(x) as.numeric(levels(x$predictions))[x$predictions])

    # Add sites with NAs
    predictions_full <- matrix(NA, nrow = nrow(vars_full), ncol = length(ranger_models))
    colnames(predictions_full) <- colnames(predictions)
    predictions_full[-inds_na,] <- predictions
    """
end
@rget predictions predictions_full inds_na

## Create Y matrices
# Get matrix Y
Y = replace(predictions_full, missing => NaN)
# Set values to NaN if no species present
inds_zeros = findall(map(x -> all(iszero.(x)), eachrow(Y)))
Y[inds_zeros,:] .= NaN
# Get Yobs
inds_obs = findall(map(x -> !any(isnan.(x)), eachrow(Y)))
inds_notobs = findall(map(x -> any(isnan.(x)), eachrow(Y)))
Yobs = Y[inds_obs,:]
# Apply Hellinger transformation
@rput Yobs
begin
    R"""
        library(vegan)
        Ytransf <- decostand(Yobs, "hel")
    """
end
@rget Ytransf

## Create distributions

# Load raw distributions (for grid size)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions

# Create RF distribution layers
Ydistrib = replace(Y, 0.0 => NaN)
rf_grids = [reshape(Ydistrib[:,i], size(distributions[1].grid)) for i in 1:size(Ydistrib, 2)]
# map(x -> reshape(Ydistrib[:,x], size(distributions[1]), 1:size(Ydistrib, 2)))
distributions = SimpleSDMResponse.(rf_grids, distributions[1].left, distributions[1].right,
                                         distributions[1].bottom, distributions[1].top)

## Export results
# save_data = true
if @isdefined save_data && save_data == true
    @save joinpath("data", "jld2", "rf-distributions.jld2") distributions
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])

## Visualize results

# Distributions
dist_rf = plotSDM(distributions[1], c=:BuPu, title = "Distributions - RF", dpi = 300)
# Richness
rich_rf = plotSDM(richness, c=:viridis, title = "Richness - RF", dpi=300)
rich_qrf = plotSDM(quantiles(richness), c=:viridis, title = "Richness quantiles - RF", dpi = 300)
# LCBD
lcbdtr_rf = plotSDM(lcbd, c=:viridis, title = "lcbd dbtr - RF", dpi = 300)
lcbdtr_qrf = plotSDM(quantiles(lcbd), c=:viridis, title = "lcbd quantiles dbtr - RF", dpi = 300)

## Export RF figures
if (@isdefined save_figures) && save_figures == true
    savefig(dist_rf, joinpath("fig", "rf", "01_rf_sp-Setophaga-coronata.png"))
    savefig(rich_rf, joinpath("fig", "rf", "03_rf_richness.png"))
    savefig(rich_qrf, joinpath("fig", "quantiles", "03_rf_richness_quantiles.png"))
    savefig(lcbdtr_rf, joinpath("fig", "rf", "05_rf_lcbd_transf.png"))
    savefig(lcbdtr_qrf, joinpath("fig", "quantiles", "05_rf_lcbd_transf_quantiles.png"))
end
