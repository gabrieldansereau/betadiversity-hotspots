import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Conditional arguments
# save_relfigures = true # should relationship figures be overwritten (optional)

## Load raw LCBD & richness results
outcome = "raw"
save_figures = false
# Load richness script
@time include("03_richness.jl")
# Load LCBD script
@time include("05_lcbd.jl")

# Keep results
raw = (distributions = distributions,
       Y = Y,
       Yobs = Yobs,
       Ytransf = Ytransf,
       inds_obs = inds_obs,
       inds_notobs = inds_notobs,
       richness = richness,
       LCBD = LCBD)

## Load SDM LCBD & richness results
outcome = "sdm"
save_figures = false
# Load richness script
@time include("03_richness.jl")
# Load LCBD script
@time include("05_lcbd.jl")

# Keep results
sdm = (distributions = distributions,
       Y = Y,
       Yobs = Yobs,
       Ytransf = Ytransf,
       inds_obs = inds_obs,
       inds_notobs = inds_notobs,
       richness = richness,
       LCBD = LCBD)

## Richness-LCBD relationship
# Calculate relative richness (α/γ)
rel_richness = [res.richness.grid ./ (size(res.Y, 2)+1) for res in (raw,sdm)]
# Scatterplot LCBD ~ richness
relation_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[1].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relation_plot, vec(rel_richness[2]), vec(sdm.LCBD[1].grid),
         markersize = 2, color=:orange, msw = 0, label = "SDM predictions")
relationtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[2].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relationtr_plot, vec(rel_richness[2]), vec(sdm.LCBD[1].grid),
         markersize = 3, c = :orange, msw = 0, label = "SDM predictions")
relationdbtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[2].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relationdbtr_plot, vec(rel_richness[2]), vec(sdm.LCBD[2].grid),
         markersize = 3, c = :orange, msw = 0, label = "SDM predictions")

## Save result
# save_relfigures = true # should relationship figures be overwritten (optional)
if (@isdefined save_relfigures) && save_relfigures == true
    @info "Figures saved (relationship)"
    savefig(relation_plot, "fig/06_relationship_lcbd-richness.png")
    savefig(relationtr_plot, "fig/06_relationship_lcbd-richness-transf.png")
else
    @info "Figures not saved (relationship)"
end

## Quantile relationship plots
qfinder = ecdf(filter(!isnan, vec(rel_richness[2])))
relation_qplot = scatter(qfinder(vec(rel_richness[1])), quantiles(vec(raw.LCBD[1].grid)),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (quantile)", ylabel = "LCBD quantile",
         grid=:none)
scatter!(relation_qplot, quantiles(vec(rel_richness[2])), quantiles(vec(sdm.LCBD[1].grid)),
         markersize = 2, color=:orange, msw = 0, label = "SDM predictions")
relationtr_qplot = scatter(qfinder(vec(rel_richness[1])), quantiles(vec(raw.LCBD[2].grid)),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (quantile)", ylabel = "LCBD quantile (raw data hellinger transformed)",
         grid=:none)
scatter!(relationtr_qplot, qfinder(vec(rel_richness[2])), quantiles(vec(sdm.LCBD[1].grid)),
         markersize = 3, c = :orange, msw = 0, label = "SDM predictions")
relationdbtr_qplot = scatter(qfinder(vec(rel_richness[1])), quantiles(vec(raw.LCBD[2].grid)),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlabel = "Species richness (quantile)", ylabel = "LCBD quantile (hellinger transformed)",
         grid=:none)
scatter!(relationdbtr_qplot, qfinder(vec(rel_richness[2])), quantiles(vec(sdm.LCBD[2].grid)),
         markersize = 3, c = :orange, msw = 0, label = "SDM predictions")
relationdbtr_nbsp_plot = scatter(vec(raw.richness.grid), quantiles(vec(raw.LCBD[2].grid)),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlabel = "Species richness (total species)", ylabel = "LCBD quantile (hellinger transformed)",
         grid=:none)
scatter!(relationdbtr_nbsp_plot, vec(sdm.richness.grid), quantiles(vec(sdm.LCBD[2].grid)),
         markersize = 3, c = :orange, msw = 0, label = "SDM predictions")

## Save result
# save_relfigures = true # should relationship figures be overwritten (optional)
if (@isdefined save_relfigures) && save_relfigures == true
    @info "Figures saved (relationship)"
    savefig(relationtr_qplot, "fig/quantiles/06_relationship_lcbd-richness-transf_quantiles.png")
    savefig(relationdbtr_nbsp_plot, "fig/quantiles/06_relationship_lcbd-richness-dbtransf-nbsp_quantiles.png")
else
    @info "Figures not saved (relationship)"
end
