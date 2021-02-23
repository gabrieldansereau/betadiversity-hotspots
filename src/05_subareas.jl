if !(@isdefined BetadiversityHotspots)
    import Pkg; Pkg.activate(".")
    @time include("required.jl")
end

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

## Subarea figures
if outcome == "raw"
    plotfct = :plotSDM2
else
    plotfct = :plot
end
function plot_lcbd_relationship(richness, lcbd; maintitle = "", kw...)
    p1 = eval(plotfct)(lcbd, c = :viridis, title = "LCBD", colorbar_title = "Relative LCBD score", clim = (0,1))
    p2 = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
                xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                xlim = (1, 50), ylim = (0.0, 1.0), clim = (1, 450),
                bottommargin = 4.0mm
                )
    if maintitle != ""
        l = @layout [t{.01h}; grid(1,2)]
        ptitle = plot(annotation = (0.5, 0.5, "$maintitle"), framestyle = :none)
        p = plot(ptitle, p1, p2, layout = l, size = (900, 300); kw...)
    else
        l = @layout [a b]
        p = plot(p1, p2, layout = l, size = (900, 300); kw...)
    end
    return p
end
resNEtr = plot_lcbd_relationship(richness_NE, lcbd_NE,
            maintitle = "Northeast subarea")
resSWtr = plot_lcbd_relationship(richness_SW, lcbd_SW,
            maintitle = "Southwest subarea")

# Combine figures
combined_plot = plot(resNEtr, resSWtr, layout = grid(2,1), 
                     size = (900, 600), 
                     bottommargin = 1.0mm,
                     title = ["" "" "" ""])

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(combined_plot, dpi = 200), joinpath("fig", outcome, "05-1_$(outcome)_subareas_combined.png"))
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
    p = plot_lcbd_relationship(richness, lcbd; kw...)
end

# Initial subarea
left = -71.0; right = -64.0; bottom = 46.0; top = 50.0;
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
# Relative LCBD values
p = plot_subareas(coords_subarea, distributions; 
                  formatter = f -> "$(round(f, digits = 1))",
                  leftmargin = 4.0mm,
                  )
# Non-relative values
asp_ratio = 92.60/60.75
p = plot_subareas(coords_subarea, distributions;
                  relative = false,
                  clim = [(0.0, Inf) (-Inf, Inf)],
                  ylim = [(-Inf, Inf) (0.0, Inf)],
                  colorbar_title = ["LCBD score" "Number of sites"],
                  formatter = :plain,
                  aspect_ratio = [asp_ratio :auto]
                  )

## Expanding GIF
left = -71.0; right = -64.0; bottom = 46.0; top = 50.0;
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
                      leftmargin = 4.0mm,
                      dpi = 200)
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
mid_ind = median(1:length(subarea_plots)) |> round |> Int64
ps = subarea_plots[[1, mid_ind, end]]

# Combine 3 scales
p = plot(ps..., dpi = 200, layout = (3,1), size = (900, 900),
         title = ["LCBD" "Relationship" "" "" "" ""],
         leftmargin = 2.0mm, bottommargin = -2.0mm,
         )

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(p, joinpath("fig", outcome, "05-2_$(outcome)_subareas_3scales.png"))
end
