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
lcbd_raw_z.grid[inds_obs] = raw_zscores
lcbd_raw_z_plot = plotSDM(lcbd_raw_z, c=:viridis)
heatmap!(lcbd_raw_z_plot, clim=(0-maximum(raw_zscores), maximum(raw_zscores)))
# SSi instead of LCBD
lcbd_ssi = copy(LCBD[4])
lcbd_ssi.grid[inds_obs] = resBDpred.SSi
plotSDM(lcbd_ssi, c=:viridis)

# Plot few results
lcbd_plot1
lcbd_ecdf_viridis
lcbd_ecdf_reds
richness_plot
richness_viridis
richness_reds
=#
