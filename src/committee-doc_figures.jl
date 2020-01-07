using Distributed
using JLD2

@time @everywhere include("src/required.jl")

## Get & prepare data
@time begin
    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
end

## Get environmental data data
# WorldClim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
# Landcover data
@time lc_vars = load_landcover(lon_range, lat_range)
# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

# Plot wcvars1
wc_plot = plotSDM(wc_vars[1], c=:auto)
heatmap!(wc_plot, clim=(-10.0, 30.0), colorbar_title="Annual Mean Temperature (°C)", dpi=300)
# Plot lcvars8 (urban)
lc_plot = plotSDM(lc_vars[2], c=:auto)
heatmap!(lc_plot, colorbar_title="Crops land cover (%)", dpi=300)

## Load data
@load "data/jld2/raw-pres-abs.jld2" pres_abs
@load "data/jld2/sdm-predictions-landcover.jld2" predictions

## Plot Single species
replace!(pres_abs[13].grid, 0.0 => NaN)
singlesp_raw = plotSDM(pres_abs[13])
title!(singlesp_raw, "")
heatmap!(singlesp_raw, clim=(0.0, 1.0), colorbar_title="Probability of occurrence", dpi=300)
heatmap!(colorbar=:none)
scatter!(singlesp_raw, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=4)
singlesp_sdm = plotSDM(predictions[13])
title!(singlesp_sdm, "")
heatmap!(singlesp_sdm, clim=(0.0, 1.0), colorbar_title="Probability of occurrence", dpi=300)

# Combine figures
singlesp_plots = plot(singlesp_raw, singlesp_sdm)
                      #layout=(2,1), size=(300,480))
## Raw
# Richness
@time include("raw/03_raw_richness.jl")
richness_raw = richness_plot
title!(richness_raw, "")
heatmap!(richness_raw, clim=(0.0, 60.0), colorbar_title = "Number of species per site", dpi=300)

# LCBD
@time include("raw/05_raw_lcbd.jl")
lcbd_raw = lcbd_plot2
title!(lcbd_raw, "")
heatmap!(lcbd_raw, colorbar_title = "LCBD value (relative to maximum)", dpi=300)

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
title!(richness_sdm, "")
heatmap!(richness_sdm, clim=(0.0, 60.0), colorbar_title = "Number of species per site", dpi=300)

# LCBD
@time include("sdm/05_sdm_lcbd.jl")
lcbd_sdm = lcbd_plot1
title!(lcbd_sdm, "")
heatmap!(lcbd_sdm, colorbar_title = "LCBD value (relative to maximum)", dpi=300)

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
        label = "Raw occurrence data",
        legend = :topright,
        xlims = (0.0, 1.0), ylims = (0.0, 1.0),
        yticks = 0.0:0.20:1.0,
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)
scatter!(relation_oneplot, rel_richness_sdm, rel_lcbd_sdm,
         markersize = 1, color=:orange, label = "SDM predictions")
# plot!(relation_oneplot, aspectratio=1)

relation_oneplot_empty = scatter(
        markersize = 1,
        c = :skyblue,
        label = "Raw occurrence data",
        legend = :bottomright,
        xlims = (0.0, 1.0), ylims = (0.0, 1.0),
        yticks = 0.0:0.20:1.0,
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)

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

# env vars
savefig(wc_plot, "doc/fig/wc_temp.png")
savefig(lc_plot, "doc/fig/lc_temp.png")

# Raw
savefig(singlesp_raw, "doc/fig/01_raw_singlesp.png")
savefig(richness_raw, "doc/fig/03_raw_richness.png")
savefig(lcbd_raw, "doc/fig/05_raw_lcbd-transf.png")
savefig(relation_raw, "doc/fig/06_raw_relation.png")

# SDM
savefig(singlesp_sdm, "doc/fig/01_sdm_singlesp.png")
savefig(richness_sdm, "doc/fig/03_sdm_richness.png")
savefig(lcbd_sdm, "doc/fig/05_sdm_lcbd.png")
savefig(relation_sdm, "doc/fig/06_sdm_relation.png")

# Combined
savefig(richness_plots, "doc/fig/03_cmb_richness.png")
savefig(lcbd_plots, "doc/fig/05_cmb_lcbd.png")
savefig(relation_plots, "doc/fig/06_cmb_relation.png")
savefig(relation_oneplot, "doc/fig/06_cmb_relation-oneplot.png")
savefig(relation_oneplot_empty, "doc/fig/06_cmb_relation-oneplot-empty.png")

## Single species SDM with threshold (run if needed, needs to load whole dataset, takes a while to run)
#=
# Option 1: sdm_single-species script, fast but not exactly the same
@time include("sdm/sdm_single-species.jl")
singlesp_sdm_threshold = maps[2]
title!(singlesp_sdm_threshold, "")
heatmap!(singlesp_sdm_threshold, clim=(0.0, 1.0), colorbar_title="Probability of occurrence", dpi=300)
savefig(singlesp_sdm_threshold, "doc/fig/01_sdm_singlesp-threshold.png")

# Option 2: sdm_predictions script, exact same but soooo long since not parallelized
## Exact same prediction, but soooooo long since not parallelized
## Get & prepare data
@time begin
    # Load data from CSV files
    df = CSV.read("data/proc/ebd_warblers_prep.csv", header=true, delim="\t")
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
    # Observed coordinates range
    lon_range_obs = extrema(df.longitude)
    lat_range_obs = extrema(df.latitude)
end

## Get environmental data
# WorldClim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
# WorldClim data with different training resolutions
#=
@time wc_vars_pred = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);
@time wc_vars_train = pmap(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], 1:19);
=#
# Landcover data
@time lc_vars = load_landcover(lon_range, lat_range)
# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

## Make prediction for Yellow Warbler
@time singlesp_pred = species_bclim(warblers_occ[13], env_vars, with_threshold=true)
singlesp_sdm_threshold = plotSDM(singlesp_pred)
title!(singlesp_sdm_threshold, "")
heatmap!(singlesp_sdm_threshold, clim=(0.0, 1.0), colorbar_title="Probability of occurrence", dpi=300)
savefig(singlesp_sdm_threshold, "doc/fig/01_sdm_singlesp-threshold.png")
=#
