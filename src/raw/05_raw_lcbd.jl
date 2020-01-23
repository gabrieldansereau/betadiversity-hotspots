import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

outcome = "raw"

## Load presence-absence data for all species
@load "data/jld2/$(outcome)-pres-abs.jld2" distributions
## Load matrix Y
@load "data/jld2/$(outcome)-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

## Compute beta diversity statistics
# Load functions
include("../lib/beta-div.jl")
# Compute BD statistics on $(outcome) data
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
t_lcbd = [fill(NaN, size(distributions[1])) for LCBDi in LCBDsets]
# Fill in grid for resBDpres & resBDperm
t_lcbd[1][inds_pred] = LCBDsets[1]
t_lcbd[2][inds_pred] = LCBDsets[2]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(t_lcbd, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], c=:viridis)
title!(lcbd_plot1, "LCBD values per site (relative to maximum, $(outcome) data)")
lcbd_plot2 = plotSDM(LCBD[2], c=:viridis)
title!(lcbd_plot2, "LCBD values per site (relative to maximum, hellinger transformed)")

## Save result
#=
savefig(lcbd_plot1, "fig/$(outcome)/05_$(outcome)_lcbd.pdf")
savefig(lcbd_plot2, "fig/$(outcome)/05_$(outcome)_lcbd-transf.pdf")
=#
