import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

outcome = "sdm"

## Load LCBD & richness scripts (if not already loaded)
# Load richness script
@time include("../$(outcome)/03_$(outcome)_richness.jl")
# Load LCBD script
@time include("../$(outcome)/05_$(outcome)_lcbd.jl")

## Investigate different ways to scale LCBD results

# LCBD values in 10 quantiles
n = 10
lcbd = copy(LCBD[1])
lcbd_nonan = filter(!isnan, lcbd.grid)
q = quantile(lcbd_nonan, 0.0:1/n:1.0)
for i in length(q)-1:-1:1
    replace!(x -> q[i] <= x <= q[i+1] ? i/10 : x, lcbd.grid)
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
lcbd_raw = copy(LCBD[3])
lcbd_raw_plot = plotSDM(lcbd_raw, c=:viridis)
# LCBD zscores
lcbd_raw_z = copy(LCBD[3])
lcbd_raw_nonan = filter(!isnan, lcbd_raw_z.grid)
raw_zscores = zscore(lcbd_raw_nonan)
lcbd_raw_z.grid[inds_obs] = raw_zscores
lcbd_raw_z_plot = plotSDM(lcbd_raw_z, c=:viridis)
heatmap!(lcbd_raw_z_plot, clim=(0-maximum(raw_zscores), maximum(raw_zscores)))
# SSi instead of LCBD
lcbd_ssi = copy(LCBD[4])
lcbd_ssi.grid[inds_obs] = resBDobs.SSi
plotSDM(lcbd_ssi, c=:viridis)

# Plot few results
lcbd_plot1
lcbd_ecdf_viridis
lcbd_ecdf_reds
richness_plot
richness_viridis
richness_reds
