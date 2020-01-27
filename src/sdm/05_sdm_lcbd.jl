import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

outcome = "sdm"

## Load distributions for all species
@load "data/jld2/$(outcome)-distributions.jld2" distributions
## Load matrix Y
@load "data/jld2/$(outcome)-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

## Compute beta diversity statistics
# Load functions
include("../lib/beta-div.jl")
# Compute BD statistics on $(outcome) data
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
t_lcbd = [fill(NaN, size(distributions[1])) for LCBDi in LCBDsets]
# Fill in grids
[t_lcbd[i][inds_obs] = LCBDsets[i] for i in 1:length(t_lcbd)]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(t_lcbd, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], c=:viridis)
title!(lcbd_plot1, "$(titlecase(outcome)) LCBD values per site (relative to maximum)")
lcbd_plot2 = plotSDM(LCBD[2], c=:viridis)
title!(lcbd_plot2, "$(titlecase(outcome)) LCBD values per site (relative to maximum, hellinger transformed)")

## Save result
#=
savefig(lcbd_plot1, "fig/$(outcome)/05_$(outcome)_lcbd.pdf")
savefig(lcbd_plot2, "fig/$(outcome)/05_$(outcome)_lcbd-transf.pdf")
=#
