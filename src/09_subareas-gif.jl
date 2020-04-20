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

#### Repeat for different subareas
function plot_subareas(coords, initial_distributions; display_coords = coords, transform = true, relative = true, kw...)
  distributions = [d[coords] for d in initial_distributions]
  Y = calculate_Ymatrix(distributions)
  richness = calculate_richness(Y.Y, Y.inds_notobs, distributions)
  lcbd = calculate_lcbd(Y.Yobs, Y.Ytransf, Y.inds_obs, distributions; relative = relative)
  if display_coords != coords
    richness = richness[display_coords]
    lcbd = [l[display_coords] for l in lcbd]
  end
  if transform
    p = plot_lcbd_richness(richness, lcbd[2]; kw...)
  else
    p = plot_lcbd_richness(richness, lcbd[1]; kw...)
  end
end

# Initial subarea
left = -71.0; right = -64.0; bottom = 47.5; top = 50.0
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
p = plot_subareas(coords_subarea, distributions; formatter = f -> "$(round(f, digits = 1))")
p = plot_subareas(coords_subarea, distributions;
                  relative = false,
                  clim = [() () () ()],
                  ylim = [() () () ()],
                  formatter = f -> "$(round(f, digits = 1))")

## Expanding GIF
left = -71.0; right = -64.0; bottom = 47.5; top = 50.0
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
subarea_plots = []
@time while left > -145.0 && bottom > 20.0
  global left -= 1.0
  global bottom -= 0.5;
  coords_subarea = (left = left, right = right, bottom = bottom, top = top)
  p = plot_subareas(coords_subarea, distributions;
                    formatter = f -> "$(round(f, digits = 1))",
                    dpi = 150)
  push!(subarea_plots, p)
end

# Create GIF
anim = @animate for p in subarea_plots
    plot(p)
end
gif(anim, fps = 3)
gif(anim, joinpath("fig", outcome, "09_subareas.gif"), fps = 3)

## Focused GIF
left = -71.0; right = -64.0; bottom = 47.5; top = 50.0
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
display_coords = coords_subarea
subarea_plots = []
@time while left > -145.0 && bottom > 20.0
  global left -= 1.0
  global bottom -= 0.5;
  coords_subarea = (left = left, right = right, bottom = bottom, top = top)
  p = plot_subareas(coords_subarea, distributions;
                    display_coords = display_coords,
                    formatter = f -> "$(round(f, digits = 1))",
                    clim = [(0,25) (0,1) (0,45) (0,1)],
                    dpi = 150)
  push!(subarea_plots, p)
end

# Create GIF
anim = @animate for p in subarea_plots
    plot(p)
end
gif(anim, fps = 7)
gif(anim, joinpath("fig", outcome, "09_subareas-focused.gif"), fps = 7)


## 3 scales comparison

ntimes = 55.0

left = -71.0; right = -64.0; bottom = 47.5; top = 50.0
coords1 = (left = left , right = right, bottom = bottom, top = top)
coords2 = (left = left - median(1:ntimes), right = right, bottom = bottom - 0.5*median(1:ntimes), top = top)
coords3 = (left = left - ntimes, right = right, bottom = bottom - 0.5*ntimes, top = top)

subarea_plots = [plot_subareas(c, distributions; formatter = f -> "$(round(f, digits = 1))") for c in (coords1, coords2, coords3)]
subarea_plots = [plot_subareas(c, distributions;
                               relative = false,
                               clim = [() (0.0,Inf) () ()],
                               ylim = [() () (0.0,Inf) ()],
                               colorbar_title = ["Number of sites" "LCBD score" "Number of sites" "LCBD score" "Number of sites" "LCBD score"],
                               formatter = :plain
                               )
                               for c in (coords1, coords2, coords3)]

ps = []
p_comb = []
for p in subarea_plots
  p_lcbd = p[2][1][:plot_object]
  p_rel = p[3][1][:plot_object]
  push!(ps, p_lcbd, p_rel)
  push!(p_comb, plot(p_lcbd, p_rel, layout = (2,1)))
end

l1 = @layout [a{0.6w} b;
              c{0.6w} d;
              e{0.6w} f]
p = plot(ps..., layout = l1, size = (1000,800))

savefig(p, joinpath("fig/", outcome, "09_$(outcome)_subareas_3scales.png"))
savefig(p, joinpath("fig/", outcome, "09_$(outcome)_subareas_3scales-abs.png"))
