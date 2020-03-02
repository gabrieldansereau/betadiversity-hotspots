import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load data
spe = CSV.read("data/proc/distributions_spe.csv", header=true, delim="\t")
spa = CSV.read("data/proc/distributions_spa.csv", header=true, delim="\t")
env = CSV.read("data/proc/distributions_env.csv", header=true, delim="\t")

## Perform RandomForests
using RCall
@rput spe spa env
begin
    R"""
    library(randomForest)

    ## Predict distributions for full range
    env_full <- read.csv("data/proc/distributions_env_full.csv", header = TRUE, sep = "\t")
    spa_full <- read.csv("data/proc/distributions_spa_full.csv", header = TRUE, sep = "\t")
    vars_full <- cbind(env_full, spa_full)
    head(vars_full)

    # Load model
    load("data/proc/rf_models.RData")

    # Make predictions
    system.time(rf_pred <- sapply(rf_models, function(x) predict(x, vars_full, type = "class")))
    """
end
@rget rf_pred

Yrf = replace(rf_pred, missing => NaN,
                           "0" => 0,
                           "1" => 1)
Yrf = Array{Float64}(Yrf)

## Load distributions for all species
@load "data/jld2/raw-distributions.jld2" distributions spenames speindex
raw = (distributions = distributions,
       spenames = spenames,
       speindex = speindex)
@load "data/jld2/sdm-distributions.jld2" distributions spenames speindex
sdm = (distributions = distributions,
       spenames = spenames,
       speindex = speindex)

rf_grids = [reshape(rf_pred[:,i], size(sdm.distributions[1])) for i in 1:size(rf_pred,2)]
rf_distributions = SimpleSDMResponse.(rf_grids, sdm.distributions[1].left, sdm.distributions[1].right,
                                         sdm.distributions[1].bottom, sdm.distributions[1].top)

plotSDM(rf_distributions[24], c=:BuPu)
plotSDM(raw.distributions[24], c=:BuPu)
plotSDM(sdm.distributions[24], c=:BuPu)
