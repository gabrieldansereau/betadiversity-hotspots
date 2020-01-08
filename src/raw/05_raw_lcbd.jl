Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")
## Load presence-absence data for all species
@load "data/jld2/raw-pres-abs.jld2" pres_abs
## Load matrix Y
@load "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

## Compute beta diversity statistics
# Load functions
include("../lib/beta-div.jl")
# Compute BD statistics on raw data
resBDpres = BD(Ypred)
# Compute BD statistics on transformed data
resBDtransf = BD(Ytransf)

# Extract LCBD values
resBD = [resBDpres, resBDtransf]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(pres_abs[1])) for LCBDi in LCBDsets]
# Fill in grid for resBDpres & resBDperm
t_lcbd[1][inds_pred] = LCBDsets[1]
t_lcbd[2][inds_pred] = LCBDsets[2]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(t_lcbd, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], c=:viridis)
title!(lcbd_plot1, "LCBD values per site (relative to maximum, raw data)")
lcbd_plot2 = plotSDM(LCBD[2], c=:viridis)
title!(lcbd_plot2, "LCBD values per site (relative to maximum, hellinger transformed)")

## Save result
#=
savefig(lcbd_plot1, "fig/raw/05_raw_lcbd.pdf")
savefig(lcbd_plot2, "fig/raw/05_raw_lcbd-transf.pdf")
=#
