using Distributed
using JLD2
using Random
@time include("required.jl")

## Load presence-absence data for all species
@load "../data/pres-abs-ebd.jld2" pres_abs
## Load matrix Y
@load "../data/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
# Compute BD statistics on raw data
resBDpres = BD(Ypred)
# Compute BD statistics on transformed data
resBDtransf = BD(Ytransf)
#=
# Run permutation tests on transformed data
@time resBDperm = BDperm(Ytransf, nperm = 999, distributed = false) # 1600 sec/27 min for 999 permutations
# Export permutation results
@save "../data/raw-resBDperm-transf.jld2" resBDperm
=#
# Load permutation results
@load "../data/raw-resBDperm-transf.jld2" resBDperm

# Extract LCBD values
resBD = [resBDpres, resBDtransf, resBDperm]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(pres_abs[1])) for LCBDi in LCBDsets]
# Fill in grid for resBDpres & resBDperm
t_lcbd[1][inds_pred] = LCBDsets[1]
t_lcbd[2][inds_pred] = LCBDsets[2]
# Get indices of sites with significant LCBDs
inds_signif = falses(size(Y,1))
inds_signif[inds_pred] = Array{Bool}(resBD[3].pLCBD .<= 0.05)
# Get indices of sites with significant LCBDs & predictions
inds_signifpred = intersect(findall(inds_signif), inds_pred)
# Fill in grid for resBDpred
t_lcbd[3][inds_pred] .= 0.0
t_lcbd[3][inds_signifpred[3]] .= 1.0
# Create SDMLayer with LCBD values
LCBD = SDMLayer.(t_lcbd, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], type="lcbd")
title!(lcbd_plot1, "LCBD values per site (relative to maximum, raw data)")
lcbd_plot2 = plotSDM(LCBD[2], type="lcbd")
title!(lcbd_plot2, "LCBD values per site (relative to maximum, hellinger transformed)")
lcbd_plot3 = plotSDM(LCBD[3], type="lcbd")
title!(lcbd_plot3, "Significant LCBDs (hellinger transformed)")

## Save result
#=
savefig(lcbd_plot1, "fig/raw-lcbd.pdf")
savefig(lcbd_plot2, "fig/raw-lcbd-transf.pdf")
savefig(lcbd_plot3, "fig/raw-lcbd-signif.pdf")
=#
