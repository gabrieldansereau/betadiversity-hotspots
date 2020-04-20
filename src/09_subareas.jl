import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

# Make sure "outcome" is defined
outcome = "rf"
if !(@isdefined outcome)
  @warn "'outcome' not defined, must be either 'raw', 'sdm' or 'rf'"
elseif !(outcome in ["raw", "sdm", "rf"])
  @warn "'outcome' invalid, must be either 'raw', 'sdm' or 'rf'"
else
  @info "'outcome' currently set to '$(outcome)'"
end

## Load distribution data for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

## Extract subareas
# Northeast subarea
coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
distributions_NE = [d[coords_NE] for d in distributions]
# Southwest subarea
coords_SW = (left = -120.0, right = -100.0, bottom = 30.0, top = 40.0)
distributions_SW = [d[coords_SW] for d in distributions]

## Get Ymatrices
NE = calculate_Ymatrix(distributions_NE)
SW = calculate_Ymatrix(distributions_SW)

## Richness
richness_NE = calculate_richness(NE.Y, NE.inds_notobs, distributions_NE)
richness_SW = calculate_richness(SW.Y, SW.inds_notobs, distributions_SW)

## LCBD
# Load functions
lcbd_NE = calculate_lcbd(NE.Yobs, NE.Ytransf, NE.inds_obs, distributions_NE)
lcbd_SW = calculate_lcbd(SW.Yobs, SW.Ytransf, SW.inds_obs, distributions_SW)

## Combine figures
function plot_lcbd_richness(richness, lcbd; title = "", kw...)
  p1 = plot(richness, c = :viridis, title = "Richness", colorbar_title = "Number of species")
  p2 = plot(lcbd, c = :viridis, title = "LCBD", colorbar_title = "Relative LCBD score", clim = (0,1))
  p3 = plot(quantiles(lcbd), c = :viridis, title = "LCBD quantiles", colorbar_title = "Quantile rank", clim = (0,1))
  p4 = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
            xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
            xlim = (1, 40), ylim = (0.0, 1.0))
  if title != ""
    l = @layout [t{.01h}; grid(2,2)]
    ptitle = plot(annotation = (0.5, 0.5, "$title"), framestyle = :none)
    p = plot(ptitle, p1, p2, p4, p3, layout = l; kw...)
  else
    p = plot(p1, p2, p4, p3; kw...)
  end
  return p
end

resNE   = plot_lcbd_richness(richness_NE, lcbd_NE[1], dpi = 150,
            title = "NE subarea - $(uppercase(outcome)) results")
resNEtr = plot_lcbd_richness(richness_NE, lcbd_NE[2], dpi = 150,
            title = "NE subarea - $(uppercase(outcome)) results (hell.transf)")

resSW   = plot_lcbd_richness(richness_SW, lcbd_SW[1], dpi = 150,
            title = "SW subarea - $(uppercase(outcome)) results")
resSWtr = plot_lcbd_richness(richness_SW, lcbd_SW[2], dpi = 150,
            title = "SW subarea - $(uppercase(outcome)) results (hell.transf)")

# Export figures
#=
savefig(resNE, joinpath("fig", outcome, "09_$(outcome)_subareas_NE.png"))
savefig(resNEtr, joinpath("fig", outcome, "09_$(outcome)_subareas_NEtr.png"))
savefig(resSW, joinpath("fig", outcome, "09_$(outcome)_subareas_SW.png"))
savefig(resSWtr, joinpath("fig", outcome, "09_$(outcome)_subareas_SWtr.png"))
=#
