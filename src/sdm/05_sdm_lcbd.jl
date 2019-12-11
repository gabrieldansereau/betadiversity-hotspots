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
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(predictions[1])) for LCBDi in LCBDsets]
# Fill in grid for resBDpred & resBDperm
t_lcbd[1][inds_pred] = LCBDsets[1]
t_lcbd[2][inds_pred] = LCBDsets[2]
t_lcbd[3][inds_pred] = LCBDsets[3]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], type="lcbd")
title!(lcbd_plot1, "SDM LCBD values per site (relative to maximum)")
lcbd_plot2 = plotSDM(LCBD[2], type="lcbd")
title!(lcbd_plot2, "SDM LCBD values per site (relative to maximum, hellinger transformed)")
lcbd_plot3 = plotSDM(LCBD[3], type="lcbd")
title!(lcbd_plot3, "SDM LCBD values per site (relative to maximum, probability data)")

## Save result
#=
savefig(lcbd_plot1, "fig/sdm/05_sdm_lcbd.pdf")
savefig(lcbd_plot2, "fig/sdm/05_sdm_lcbd-transf.pdf")
savefig(lcbd_plot3, "fig/sdm/05_sdm_lcbd-prob.pdf")
=#
