using Distributed
using JLD2
using Random
@time include("required.jl")

## Load predictions for all species
@load "data/jld2/predictions-ebd.jld2" predictions
## Load matrix Y
@load "data/jld2/sdm-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
# Compute BD statistics on predictions
resBDpred = BD(Ypred)
# Compute BD statistics on transformed predictions
resBDtransf = BD(Ytransf)
#=
# Run permutation tests on transformed predictions
@time resBDperm = BDperm(Ytransf, nperm = 999, distributed = false) # 300 sec/5 min for 999 permutations on Ada
# Export permutation results
@save "data/jld2/sdm-resBDperm-transf.jld2" resBDperm
=#
# Load permutation results
@load "data/jld2/sdm-resBDperm-transf.jld2" resBDperm

# Extract LCBD values
resBD = [resBDpred, resBDtransf, resBDperm]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(predictions[1])) for LCBDi in LCBDsets]
# Fill in grid for resBDpred & resBDperm
t_lcbd[1][inds_pred] = LCBDsets[1]
t_lcbd[2][inds_pred] = LCBDsets[2]
# Get indices of sites with significant LCBDs
inds_signif = falses(size(Y,1))
inds_signif[inds_pred] = Array{Bool}(resBD[3].pLCBD .<= 0.05)
# Get indices of sites with significant LCBDs & predictions
inds_signifpred = intersect(findall(inds_signif), inds_pred)
# Fill in grid for resBDpred
t_lcbd[3][inds_pred] .= 0.0
t_lcbd[3][inds_signifpred] .= 1.0
# Create SDMLayer with LCBD values
LCBD = SDMLayer.(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], type="lcbd")
title!(lcbd_plot1, "SDM LCBD values per site (relative to maximum)")
lcbd_plot2 = plotSDM(LCBD[2], type="lcbd")
title!(lcbd_plot2, "SDM LCBD values per site (relative to maximum, hellinger transformed)")
lcbd_plot3 = plotSDM(LCBD[3], type="lcbd")
title!(lcbd_plot3, "SDM significant LCBDs (hellinger transformed)")

## Save result
#=
savefig(lcbd_plot1, "fig/sdm/sdm-lcbd.pdf")
savefig(lcbd_plot2, "fig/sdm/sdm-lcbd-transf.pdf")
savefig(lcbd_plot3, "fig/sdm/sdm-lcbd-signif.pdf")
=#
