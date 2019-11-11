using Distributed
using JLD2
using Random
@time include("../required.jl")

## Load predictions for all species
@load "data/jld2/sdm-predictions.jld2" predictions
## Load matrix Y
@load "data/jld2/sdm-Y-matrices.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred

## Compute beta diversity statistics
# Load functions
include("../lib/beta-div.jl")
#=
# Run permutation tests on transformed predictions
@time resBDperm = BDperm(Ytransf, nperm = 999, distributed = false) # 300 sec/5 min for 999 permutations on Ada
# Export permutation results
@save "data/jld2/sdm-resBDperm-transf.jld2" resBDperm
=#
# Load permutation results
@load "data/jld2/sdm-resBDperm-transf.jld2" resBDperm

# Extract LCBD values
resBD = [resBDperm]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(predictions[1])) for LCBDi in LCBDsets]
# Get indices of sites with significant LCBDs
inds_signif = falses(size(Y,1))
inds_signif[inds_pred] = Array{Bool}(resBD[1].pLCBD .<= 0.05)
# Get indices of sites with significant LCBDs & predictions
inds_signifpred = intersect(findall(inds_signif), inds_pred)
# Fill in grid for resBDpred
t_lcbd[1][inds_pred] .= 0.0
t_lcbd[1][inds_signifpred] .= 1.0
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
lcbd_plot4 = plotSDM(LCBD[1], type="lcbd")
title!(lcbd_plot4, "SDM significant LCBDs (hellinger transformed)")

## Save result
#=
savefig(lcbd_plot4, "fig/sdm/05_sdm_lcbd-signif.pdf")
=#