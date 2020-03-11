import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load data
#=
spe = CSV.read("data/proc/distributions_spe.csv", header=true, delim="\t")
spa = CSV.read("data/proc/distributions_spa.csv", header=true, delim="\t")
env = CSV.read("data/proc/distributions_env.csv", header=true, delim="\t")
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

# Create Y matrices
Yrf = replace(predictions_full, missing => NaN)
Yrf = Array{Float64}(Yrf)
Yobs = predictions
inds_notobs = inds_na
inds_obs = collect(1:size(Yrf,1))[Not(inds_notobs)]

# Load raw distributions (for grid size)
@load "data/jld2/raw-distributions.jld2" distributions

# Create RF distribution layers
Ydistrib = replace(Yrf, 0.0 => NaN)
rf_grids = [reshape(Ydistrib[:,i], size(distributions[1])) for i in 1:size(Ydistrib, 2)]
# map(x -> reshape(Ydistrib[:,x], size(distributions[1]), 1:size(Ydistrib, 2)))
rf_distributions = SimpleSDMResponse.(rf_grids, distributions[1].left, distributions[1].right,
                                         distributions[1].bottom, distributions[1].top)

plotSDM(rf_distributions[1], c=:BuPu)

#### Richness ####

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Yrf))
# Reshape to grid format
sums = reshape(sums, size(distributions[1]))

## Create SimpleSDMLayer
richness = SimpleSDMResponse(sums, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

#### LCBD

inds_obs = findall(map(x -> !any(isnan.(x)), eachrow(Yrf)))
inds_notobs = findall(map(x -> any(isnan.(x)), eachrow(Yrf)))
Yobs = Yrf[inds_obs,:]
@rput Yobs
begin
    R"""
        library(vegan)
        Ytransf <- decostand(Yobs, "hel")
    """
end
@rget Ytransf

# Load functions
include("lib/beta-div.jl")
# Compute BD statistics on distribution data
resBDobs = BD(Yobs)
# Compute BD statistics on transformed data
resBDtransf = BD(Ytransf)

# Extract LCBD values
resBD = [resBDobs, resBDtransf]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets_raw = copy(LCBDsets)
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]
LCBDsets = [LCBDsets..., LCBDsets_raw...]

## Arrange LCBD values as grid
# Create empty grids
LCBDgrids = [fill(NaN, size(distributions[1])) for LCBDi in LCBDsets]
# Fill in grids
[LCBDgrids[i][inds_obs] = LCBDsets[i] for i in 1:length(LCBDgrids)]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(LCBDgrids, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

#### Compare with previous results
## Save random forest results
rf = (distributions = rf_distributions,
      Y = Yrf,
      Yobs = Yobs,
      Ytransf = Ytransf,
      inds_obs = inds_obs,
      inds_notobs = inds_notobs,
      richness = richness,
      LCBD = LCBD)

## Load raw LCBD & richness results
outcome = "raw"
save_figures = false
# Load richness script
@time include("03_richness.jl")
# Load LCBD script
@time include("05_lcbd.jl")

# Keep results
raw = (distributions = distributions,
       Y = Y,
       Yobs = Yobs,
       Ytransf = Ytransf,
       inds_obs = inds_obs,
       inds_notobs = inds_notobs,
       richness = richness,
       LCBD = LCBD)

## Load SDM LCBD & richness results
outcome = "sdm"
save_figures = false
# Load richness script
@time include("03_richness.jl")
# Load LCBD script
@time include("05_lcbd.jl")

# Keep results
sdm = (distributions = distributions,
       Y = Y,
       Yobs = Yobs,
       Ytransf = Ytransf,
       inds_obs = inds_obs,
       inds_notobs = inds_notobs,
       richness = richness,
       LCBD = LCBD)

dist_rf = plotSDM(rf.distributions[1], c=:BuPu, title = "Distributions - RF")
dist_raw = plotSDM(raw.distributions[1], c=:BuPu, title = "Distributions - Raw")
dist_sdm = plotSDM(sdm.distributions[1], c=:BuPu, title = "Distributions - SDM")

rich_rf = plotSDM(rf.richness, c=:viridis, title = "Richness - RF")
rich_raw = plotSDM(raw.richness, c=:viridis, title = "Richness - Raw")
rich_sdm = plotSDM(sdm.richness, c=:viridis, title = "Richness - SDM")

rich_qrf = plotSDM(quantiles(rf.richness), c=:viridis, title = "Richness quantiles - RF")
rich_qraw = plotSDM(quantiles(raw.richness), c=:viridis, title = "Richness quantiles - Raw")
rich_qsdm = plotSDM(quantiles(sdm.richness), c=:viridis, title = "Richness quantiles - SDM")

lcbd_rf = plotSDM(rf.LCBD[1], c=:viridis, title = "LCBD - RF")
lcbd_raw = plotSDM(raw.LCBD[1], c=:viridis, title = "LCBD - Raw")
lcbd_sdm = plotSDM(sdm.LCBD[1], c=:viridis, title = "LCBD - SDM")

lcbdtr_rf = plotSDM(rf.LCBD[2], c=:viridis, title = "LCBD dbtr - RF")
lcbdtr_raw = plotSDM(raw.LCBD[2], c=:viridis, title = "LCBD dbtr - Raw")
lcbdtr_sdm = plotSDM(sdm.LCBD[2], c=:viridis, title = "LCBD dbtr - SDM")

lcbd_qrf = plotSDM(quantiles(rf.LCBD[1]), c=:viridis, title = "LCBD quantiles - RF")
lcbd_qraw = plotSDM(quantiles(raw.LCBD[1]), c=:viridis, title = "LCBD quantiles - Raw")
lcbd_qsdm = plotSDM(quantiles(sdm.LCBD[1]), c=:viridis, title = "LCBD quantiles - SDM")

lcbdtr_qrf = plotSDM(quantiles(rf.LCBD[2]), c=:viridis, title = "LCBD quantiles dbtr - RF")
lcbdtr_qraw = plotSDM(quantiles(raw.LCBD[2]), c=:viridis, title = "LCBD quantiles dbtr - Raw")
lcbdtr_qsdm = plotSDM(quantiles(sdm.LCBD[2]), c=:viridis, title = "LCBD quantiles dbtr - SDM")

dist = plot(dist_rf, dist_raw, dist_sdm)
rich = plot(rich_rf, rich_raw, rich_sdm)
rich_q = plot(rich_qrf, rich_qraw, rich_qsdm)
lcbd = plot(lcbd_rf, lcbd_raw, lcbd_sdm)
lcbdtr = plot(lcbdtr_rf, lcbdtr_raw, lcbdtr_sdm)
lcbd_q = plot(lcbd_qrf, lcbd_qraw, lcbd_qsdm)
lcbdtr_q = plot(lcbdtr_qrf, lcbdtr_qraw, lcbdtr_qsdm)

dist
rich
rich_q
lcbd
lcbdtr
lcbd_q
lcbdtr_q

#### Relationship
# Calculate relative richness (α/γ)
rel_richness = [res.richness.grid ./ (size(res.Y, 2)+1) for res in (raw,rf)]
# Scatterplot LCBD ~ richness
relation_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[1].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relation_plot, vec(rel_richness[2]), vec(rf.LCBD[1].grid),
         markersize = 2, color=:orange, msw = 0, label = "SDM predictions")
relationtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[2].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relationtr_plot, vec(rel_richness[2]), vec(rf.LCBD[1].grid),
         markersize = 3, c = :orange, msw = 0, label = "SDM predictions")
relationdbtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[2].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relationdbtr_plot, vec(rel_richness[2]), vec(rf.LCBD[2].grid),
         markersize = 2, c = :orange, msw = 0, label = "SDM predictions")
