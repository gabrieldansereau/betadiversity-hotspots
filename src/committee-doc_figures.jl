using Distributed
using JLD2

@time @everywhere include("src/required.jl")

## Get & prepare data
@time begin
    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);

## Load data
@load "data/jld2/raw-pres-abs.jld2" pres_abs
@load "data/jld2/sdm-predictions.jld2" predictions

## Plot Single species
replace!(pres_abs[13].grid, 0.0 => NaN)
singlesp_raw = plotSDM(pres_abs[13])
title!(singlesp_raw, "")
heatmap!(singlesp_raw, clim=(0.0, 1.0), colorbar_title="Probability of occurrence")
singlesp_sdm = plotSDM(predictions[13])
heatmap!(singlesp_sdm, clim=(0.0, 1.0), colorbar_title="Probability of occurrence")
title!(singlesp_sdm, "")

# Combine figures
singlesp_plots = plot(singlesp_raw, singlesp_sdm)
                      #layout=(2,1), size=(300,480))
## Raw
# Richness
@time include("raw/03_raw_richness.jl")
richness_raw = richness_plot
heatmap!(richness_raw, clim=(0.0, 60.0), colorbar_title = "Number of species per site")
title!(richness_raw, "")

# LCBD
@time include("raw/05_raw_lcbd.jl")
lcbd_raw = lcbd_plot2
title!(lcbd_raw, "")
heatmap!(lcbd_raw, colorbar_title = "LCBD value (relative to maximum)")

# Relationship
@time include("raw/06_raw_relation-lcbd-richness.jl")
relation_raw = relationtr_plot
title!(relation_raw, "")
xaxis!(relation_raw, (0.0, 1.0))
yaxis!(relation_raw, (0.0, 1.0))

# Keep richness & lcbd values for combined relationship plot
rel_richness_raw = vec(rel_richness)
rel_lcbd_raw = vec(LCBD[2].grid)

## SDM

# Richness
@time include("sdm/03_sdm_richness.jl")
richness_sdm = richness_plot
heatmap!(richness_sdm, clim=(0.0, 60.0), colorbar_title = "Number of species per site")
title!(richness_sdm, "")

# LCBD
@time include("sdm/05_sdm_lcbd.jl")
lcbd_sdm = lcbd_plot1
title!(lcbd_sdm, "")
heatmap!(lcbd_sdm, colorbar_title = "LCBD value (relative to maximum)")

# Relationship
@time include("sdm/06_sdm_relation-lcbd-richness.jl")
relation_sdm = relation_plot
title!(relation_sdm, "")
xaxis!(relation_sdm, (0.0, 1.0))
yaxis!(relation_sdm, (0.0, 1.0))

# Keep richness & lcbd values for combined relationship plot
rel_richness_sdm = vec(rel_richness)
rel_lcbd_sdm = vec(LCBD[1].grid)

## Combine plots
# Side by side
richness_plots = plot(richness_raw, richness_sdm)
lcbd_plots = plot(lcbd_raw, lcbd_sdm)
relation_plots = plot(relation_raw, relation_sdm)

# Relationships on same plot
relation_oneplot = scatter(rel_richness_raw, rel_lcbd_raw,
        markersize = 1,
        c = :skyblue,
        label = "Raw",
        legend = :bottomright,
        xlims = (0.0, 1.0), ylims = (0.0, 1.0),
        yticks = 0.0:0.20:1.0,
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)
scatter!(relation_oneplot, rel_richness_sdm, rel_lcbd_sdm,
         markersize = 1, color=:orange, label = "SDM")
# plot!(relation_oneplot, aspectratio=1)

# Mean LCBD / richness
rel_df_raw = DataFrame(richness = rel_richness_raw, lcbd = rel_lcbd_raw)
rel_mean_raw = sort(by(rel_df_raw, :richness, lcbd_mean = :lcbd => mean), :richness)
rel_mean_plot = scatter(rel_mean_raw.richness, rel_mean_raw.lcbd_mean)

rel_df_sdm = DataFrame(richness = rel_richness_sdm, lcbd = rel_lcbd_sdm)
rel_mean_sdm = sort(by(rel_df_sdm, :richness, lcbd_mean = :lcbd => mean), :richness)
scatter!(rel_mean_plot, rel_mean_sdm.richness, rel_mean_sdm.lcbd_mean)

# Median LCBD / richness
rel_df_raw = DataFrame(richness = rel_richness_raw, lcbd = rel_lcbd_raw)
rel_median_raw = sort(by(rel_df_raw, :richness, lcbd_median = :lcbd => median), :richness)
rel_median_plot = plot(rel_median_raw.richness, rel_median_raw.lcbd_median, smooth=true)

rel_df_sdm = DataFrame(richness = rel_richness_sdm, lcbd = rel_lcbd_sdm)
rel_median_sdm = sort(by(rel_df_sdm, :richness, lcbd_median = :lcbd => median), :richness)
plot!(rel_median_plot, rel_median_sdm.richness, rel_median_sdm.lcbd_median, smooth=true)

# Density
#=
using StatsPlots
density_richness_raw = density(filter(!isnan, rel_richness_raw))
density_richness_sdm = density(filter(!isnan, rel_richness_sdm))

density_lcbd_raw = density(filter(!isnan, rel_lcbd_raw))
density_lcbd_sdm = density(filter(!isnan, rel_lcbd_sdm))
=#

## Export figures
# Raw
savefig(singlesp_raw, "doc/fig/01_raw_singlesp.pdf")
savefig(richness_raw, "doc/fig/03_raw_richness.pdf")
savefig(lcbd_raw, "doc/fig/05_raw_lcbd-transf.pdf")
savefig(relation_raw, "doc/fig/06_raw_relation.png")

# Richness
savefig(singlesp_sdm, "doc/fig/01_sdm_singlesp.pdf")
savefig(richness_sdm, "doc/fig/03_sdm_richness.pdf")
savefig(lcbd_sdm, "doc/fig/05_sdm_lcbd.pdf")
savefig(relation_sdm, "doc/fig/06_sdm_relation.png")

# Combined
savefig(richness_plots, "doc/fig/03_cmb_richness.pdf")
savefig(lcbd_plots, "doc/fig/05_cmb_lcbd.pdf")
savefig(relation_plots, "doc/fig/06_cmb_relation.png")
savefig(relation_oneplot, "doc/fig/06_cmb_relation-oneplot.png")
