using Distributed
using JLD2
using Random
@time include("../required.jl")

## Load predictions for all species
@load "data/jld2/sdm-predictions-landcover.jld2" predictions
## Load matrix Y
@load "data/jld2/sdm-Y-matrices-landcover.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred

## Compute beta diversity statistics
# Load functions
include("../lib/beta-div.jl")
# Compute BD statistics on binary predictions
resBDpred = BD(Ypred)
# Compute BD statistics on transformed binary predictions
resBDtransf = BD(Ytransf)
# Compute BD statistics on probability predictions
resBDprob = BD(Yprob)

# Extract LCBD values
resBD = [resBDpred, resBDtransf, resBDprob]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets_raw = copy(LCBDsets)
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]
LCBDsets = [LCBDsets..., LCBDsets_raw...]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(predictions[1])) for LCBDi in LCBDsets]
# Fill in grid for resBDpred & resBDperm
[t_lcbd[i][inds_pred] = LCBDsets[i] for i in 1:length(t_lcbd)]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], c=:viridis)
title!(lcbd_plot1, "SDM LCBD values per site (relative to maximum)")
lcbd_plot2 = plotSDM(LCBD[2], c=:viridis)
title!(lcbd_plot2, "SDM LCBD values per site (relative to maximum, hellinger transformed)")
lcbd_plot3 = plotSDM(LCBD[3], c=:viridis)
title!(lcbd_plot3, "SDM LCBD values per site (relative to maximum, probability data)")

## Save result
#=
savefig(lcbd_plot1, "fig/sdm/05_sdm_lcbd.pdf")
savefig(lcbd_plot2, "fig/sdm/05_sdm_lcbd-transf.pdf")
savefig(lcbd_plot3, "fig/sdm/05_sdm_lcbd-prob.pdf")
=#


## Investigate different ways to scale LCBD results

#=
# LCBD values in 10 quantiles
n = 10
lcbd = copy(LCBD[1])
lcbd_nonan = filter(!isnan, lcbd.grid)
q = quantile(lcbd_nonan, 0.0:1/n:1.0)
for i in length(q)-1:-1:1
    replace!(x -> q[i] <= x <= q[i+1] ? i : x, lcbd.grid)
end
by(DataFrame(lcbd = vec(lcbd.grid)), :lcbd, nrow)
plotSDM(lcbd, c=:viridis)
plotSDM(lcbd, c=:magma)

# LCBD values by ecdf distribution
lcbd_ecdf = copy(LCBD[1])
qfinder = ecdf(filter(!isnan, lcbd_ecdf.grid))
replace!(x -> !isnan(x) ? qfinder(x) : x, lcbd_ecdf.grid)
lcbd_ecdf_viridis = plotSDM(lcbd_ecdf, c=:viridis)
lcbd_ecdf_reds = plotSDM(lcbd_ecdf, c=:reds)


# Raw LCBD scores
lcbd_raw = copy(LCBD[4])
lcbd_raw_plot = plotSDM(lcbd_raw, c=:viridis)
lcbd_raw_z = copy(LCBD[4])
lcbd_raw_nonan = filter(!isnan, lcbd_std.grid)
# LCBD zscores
raw_zscores = zscore(lcbd_raw_nonan)
lcbd_raw_z.grid[inds_pred] = raw_zscores
lcbd_raw_z_plot = plotSDM(lcbd_raw_z, c=:viridis)
heatmap!(lcbd_raw_z_plot, clim=(0-maximum(raw_zscores), maximum(raw_zscores)))
# SSi instead of LCBD
lcbd_ssi = copy(LCBD[4])
lcbd_ssi.grid[inds_pred] = resBDpred.SSi
plotSDM(lcbd_ssi, c=:viridis)

# Plot few results
lcbd_plot1
lcbd_ecdf_viridis
lcbd_ecdf_reds
richness_plot
richness_viridis
richness_reds
=#
