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

# Load raw distributions (for grid size)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions

# Create RF distribution layers
Ydistrib = replace(Y, 0.0 => NaN)
rf_grids = [reshape(Ydistrib[:,i], size(distributions[1].grid)) for i in 1:size(Ydistrib, 2)]
# map(x -> reshape(Ydistrib[:,x], size(distributions[1]), 1:size(Ydistrib, 2)))
distributions = SimpleSDMResponse.(rf_grids, distributions[1].left, distributions[1].right,
                                         distributions[1].bottom, distributions[1].top)

plotSDM(distributions[1], c=:BuPu)

## Richness

richness = calculate_richness(Y, inds_notobs, distributions[1])
plotSDM(richness, c = :viridis)

## LCBD

LCBD = calculate_lcbd(Yobs, Ytransf, inds_obs, distributions[1])
plotSDM(LCBD, c = :viridis)

## Export results
@save joinpath("data", "jld2", "rf-distributions.jld2") distributions
@save joinpath("data", "jld2", "rf-Y-matrices.jld2") Y Yobs Ytransf inds_obs inds_notobs
# Load results
@load joinpath("data", "jld2", "rf-distributions.jld2") distributions
@load joinpath("data", "jld2", "rf-Y-matrices.jld2") Y Yobs Ytransf inds_obs inds_notobs

#### Compare with previous results
## Save random forest results
rf = (distributions = distributions,
      Y = Y,
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

## Compare results

dist_rf = plotSDM(rf.distributions[1], c=:BuPu, title = "Distributions - RF", dpi = 300)
dist_raw = plotSDM(raw.distributions[1], c=:BuPu, title = "Distributions - Raw", dpi = 300)
dist_sdm = plotSDM(sdm.distributions[1], c=:BuPu, title = "Distributions - SDM", dpi = 300)

rich_rf = plotSDM(rf.richness, c=:viridis, title = "Richness - RF", dpi=300)
rich_raw = plotSDM(raw.richness, c=:viridis, title = "Richness - Raw", dpi = 300)
rich_sdm = plotSDM(sdm.richness, c=:viridis, title = "Richness - SDM", dpi = 300)

rich_qrf = plotSDM(quantiles(rf.richness), c=:viridis, title = "Richness quantiles - RF", dpi = 300)
rich_qraw = plotSDM(quantiles(raw.richness), c=:viridis, title = "Richness quantiles - Raw", dpi = 300)
rich_qsdm = plotSDM(quantiles(sdm.richness), c=:viridis, title = "Richness quantiles - SDM", dpi = 300)

lcbdtr_rf = plotSDM(rf.LCBD, c=:viridis, title = "LCBD dbtr - RF", dpi = 300)
lcbdtr_raw = plotSDM(raw.LCBD[2], c=:viridis, title = "LCBD dbtr - Raw", dpi = 300)
lcbdtr_sdm = plotSDM(sdm.LCBD[2], c=:viridis, title = "LCBD dbtr - SDM", dpi = 300)

lcbdtr_qrf = plotSDM(quantiles(rf.LCBD), c=:viridis, title = "LCBD quantiles dbtr - RF", dpi = 300)
lcbdtr_qraw = plotSDM(quantiles(raw.LCBD[2]), c=:viridis, title = "LCBD quantiles dbtr - Raw", dpi = 300)
lcbdtr_qsdm = plotSDM(quantiles(sdm.LCBD[2]), c=:viridis, title = "LCBD quantiles dbtr - SDM", dpi = 300)

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
relationdbtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[2].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid = :none,
         dpi = 300)
scatter!(relationdbtr_plot, vec(rel_richness[2]), vec(rf.LCBD.grid),
         markersize = 2, c = :orange, msw = 0, label = "SDM predictions")

#### Export RF figures
if (@isdefined save_figures) && save_figures == true
    savefig(dist_rf, joinpath("fig", "rf", "01_rf_sp-Setophaga-coronata.png"))
    savefig(rich_rf, joinpath("fig", "rf", "03_rf_richness.png"))
    savefig(rich_qrf, joinpath("fig", "quantiles", "03_rf_richness_quantiles.png"))
    savefig(lcbdtr_rf, joinpath("fig", "rf", "05_rf_lcbd_transf.png"))
    savefig(lcbdtr_qrf, joinpath("fig", "quantiles", "05_rf_lcbd_transf_quantiles.png"))
    savefig(relationdbtr_plot, joinpath("fig", "rf", "06_rf_relationship_dbtransf.png"))
end
