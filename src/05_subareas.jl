import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

## Conditional arguments
# outcome = "rf"
outcome = "bart"
# save_figures = true

# Make sure "outcome" is defined
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw', 'bio', 'rf', or 'bart'"
elseif !(outcome in ["raw", "bio", "rf", "bart"])
    @warn "'outcome' invalid, must be either 'raw', 'bio', 'rf', or 'bart'"
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
Y_NE = calculate_Y(distributions_NE)
Y_SW = calculate_Y(distributions_SW)

## Richness
richness_NE = calculate_richness(Y_NE, distributions_NE[1])
richness_SW = calculate_richness(Y_SW, distributions_SW[1])

## LCBD
# Load functions
lcbd_NE = calculate_lcbd(Y_NE, distributions_NE[1])
lcbd_SW = calculate_lcbd(Y_SW, distributions_SW[1])

## Combine figures
function plot_lcbd_richness(richness, lcbd; title = "", kw...)
    p1 = plot(richness, c = :viridis, title = "Richness", colorbar_title = "Number of species")
    p2 = plot(lcbd, c = :viridis, title = "LCBD", colorbar_title = "Relative LCBD score", clim = (0,1))
    p3 = plot(quantiles(lcbd), c = :viridis, title = "LCBD quantiles", colorbar_title = "Quantile rank", clim = (0,1))
    p4 = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
                     xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                     xlim = (1, 45), ylim = (0.0, 1.0))
    if title != ""
        l = @layout [t{.01h}; grid(2,2)]
        ptitle = plot(annotation = (0.5, 0.5, "$title"), framestyle = :none)
        p = plot(ptitle, p1, p2, p4, p3, layout = l; kw...)
    else
        p = plot(p1, p2, p4, p3; kw...)
    end
    return p
end

resNEtr = plot_lcbd_richness(richness_NE, lcbd_NE,
            title = "NE subarea - $(uppercase(outcome)) results (hell.transf)")
resSWtr = plot_lcbd_richness(richness_SW, lcbd_SW,
            title = "SW subarea - $(uppercase(outcome)) results (hell.transf)")

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(resNEtr, dpi = 150), joinpath("fig", outcome, "05-1_$(outcome)_subareas_NEtr.png"))
    savefig(plot(resSWtr, dpi = 150), joinpath("fig", outcome, "05-1_$(outcome)_subareas_SWtr.png"))
end

#### Repeat for different subareas
function plot_subareas(coords, initial_distributions; display_coords = coords, transform = true, relative = true, kw...)
    distributions = [d[coords] for d in initial_distributions]
    Y = calculate_Y(distributions)
    richness = calculate_richness(Y, distributions[1])
    lcbd = calculate_lcbd(Y, distributions[1];
                          transform = transform, relative = relative)
    if display_coords != coords
        richness = richness[display_coords]
        lcbd = [l[display_coords] for l in lcbd]
    end
    p = plot_lcbd_richness(richness, lcbd; kw...)
end

# Initial subarea
left = -71.0; right = -64.0; bottom = 46.5; top = 50.0;
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
# Relative LCBD values
p = plot_subareas(coords_subarea, distributions; formatter = f -> "$(round(f, digits = 1))")
# Non-relative values
p = plot_subareas(coords_subarea, distributions;
                  relative = false,
                  clim = [() (0.0,Inf) () ()],
                  ylim = [() () (0.0,Inf) ()],
                  colorbar_title = ["Number of sites" "LCBD score" "Number of sites" "LCBD score" "Number of sites" "LCBD score"],
                  formatter = :plain,
                  aspect_ratio = [asp_ratio asp_ratio :auto asp_ratio]
                  )

## Expanding GIF
left = -71.0; right = -64.0; bottom = 46.5; top = 50.0;
dim_ratio = (top-bottom)/(right-left)
asp_ratio = 92.60/60.75
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
subarea_plots = []
nplots = 0
@time while left > -145.0 + asp_ratio && bottom > 20.0 + asp_ratio * dim_ratio
    global nplots += 1
    global left -= asp_ratio
    global bottom -= asp_ratio * dim_ratio
    coords_subarea = (left = left, right = right, bottom = bottom, top = top)
    p = plot_subareas(coords_subarea, distributions;
                      formatter = f -> "$(round(f, digits = 1))",
                      aspect_ratio = [asp_ratio asp_ratio :auto asp_ratio],
                      dpi = 150)
    push!(subarea_plots, p)
end

# Create GIF
anim = @animate for p in subarea_plots[Not(1)]
    plot(p)
end
gif(anim, fps = 7)
if (@isdefined save_figures) && save_figures == true
    gif(anim, joinpath("fig", outcome, "05-3_$(outcome)_subareas.gif"), fps = 3)
end

#### 3 scales comparison

# Extract LCBD & relationship subplots for first, middle, last GIF plots
ps = []
mid_ind = median(1:length(subarea_plots)) |> round |> Int64
for p in subarea_plots[[1, mid_ind, end]]
    p_lcbd = p[2][1][:plot_object]
    p_rel = p[3][1][:plot_object]
    push!(ps, p_lcbd, p_rel)
end

# Combine 3 scales
l1 = @layout [a{0.6w} b;
              c{0.6w} d;
              e{0.6w} f]
p = plot(ps..., layout = l1, size = (1000,800))

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(p, joinpath("fig", outcome, "05-2_$(outcome)_subareas_3scales.png"))
end
