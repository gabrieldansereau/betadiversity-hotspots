include("required.jl")

## Conditional arguments
# outcome = "raw"
# outcome = "bart"
# save_figures = true

# Make sure "outcome" is defined
if !(@isdefined outcome) || !(outcome in ["raw", "bart"])
    @error "'outcome' must be either 'raw' or 'bart'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Load distribution data for all species

# Load species names
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
spenames = filter(:type => ==("species"), glossary).full_name

# Load distributions
distributions = [
    geotiff(
        SimpleSDMPredictor, joinpath("data", "raster", "distributions_$(outcome).tif"), i
    ) for i in eachindex(spenames)
]

## Extract subareas

# Northeast subarea
coords_NE = (left=-80.0, right=-60.0, bottom=40.0, top=50.0)
distributions_NE = [d[coords_NE] for d in distributions]

# Southwest subarea
coords_SW = (left=-120.0, right=-100.0, bottom=30.0, top=40.0)
distributions_SW = [d[coords_SW] for d in distributions]

# Get Ymatrices
Y_NE = Ymatrix(distributions_NE)
Y_SW = Ymatrix(distributions_SW)

# Get Y matrices for observed sites only
Yobs_NE = Ymatrix(distributions_NE; observed=true)
Yobs_SW = Ymatrix(distributions_SW; observed=true)

# Richness
richness_NE = richness(Y_NE, distributions_NE[1])
richness_SW = richness(Y_SW, distributions_SW[1])

# Absolute LCBD values
lcbd_NE = lcbd(Y_NE, distributions_NE[1]; relative=false)
lcbd_SW = lcbd(Y_SW, distributions_SW[1]; relative=false)

# Total beta
beta_NE = beta_total(Y_NE)
beta_SW = beta_total(Y_SW)

## Subarea figures

# Choose plotting function according to outcome
if outcome == "raw"
    plotfct = :plot_layer # slower, but background is nice for sparse data
else
    plotfct = :plot # faster, sufficient for continuous data
end

# Functions to add total beta diversity in rectangle box
rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
function rectangle!(w, h, x, y, textstring, textsize)
    plot!(
        rectangle(w, h, x, y),
        color=:white, legend=:none,
        annotations=(x + (w / 2), y + (h / 2), text(textstring, textsize, :center)),
    )
end

# Function to produce combined plots
function plot_lcbd_relationship(
    richness_layer, lcbd_layer, beta_total; scale=true, scaling_value=1, kw...
)
    # Scale values for a nicer result
    if scale
        # Set scaling value automatically (optional, set scaling_value=... otherwise)
        if scaling_value == 1
            scaling_factor = lcbd_layer |> maximum |> log10 |> abs |> ceil |> Int
            scaling_value = 10^scaling_factor
        end
        lcbd_layer = rescale(lcbd_layer, extrema(lcbd_layer) .* scaling_value)
    end

    # LCBD map subpanel
    p1 = eval(plotfct)(
        lcbd_layer;
        c=:viridis,
        colorbar_title="LCBD value (x $(format(scaling_value, commas=true)))",
        clim=extrema(lcbd_layer)
    )

    # Relationship subpanel
    p2 = histogram2d(
        richness_layer,
        lcbd_layer;
        c=:viridis,
        bins=40,
        colorbar_title="Number of sites",
        xlabel="Richness",
        ylabel="LCBD value (x $(format(scaling_value, commas=true)))",
        xlim=(1, 50),
        ylim=(-Inf, Inf),
    )
    vline!([median(richness_layer)], label=:none, linestyle=:dash, c=:grey)
    hline!([median(lcbd_layer)], label=:none, linestyle=:dash, c=:grey)

    # Add total beta value in rectangle box (fitted to dimensions)
    lmin, lmax = extrema(lcbd_layer)
    lrange = lmax - lmin
    rectangle!(
        16.0,
        0.15 * lrange,
        33.0,
        lmax - 0.2 * lrange,
        "BDtot = $(round(beta_total; digits=3))",
        7,
    )

    # Combine in subpanels in single plot
    l = @layout [a{0.5w} b{0.5w,0.9h}]
    p = plot(
        p1,
        p2;
        layout=l,
        size=(900, 300),
        # bottommargin = [0.0mm 5.0mm],
        rightmargin=[5.0mm 0mm],
        leftmargin=[5.0mm 5.0mm],
        kw...
    )
    return p
end

# Create subarea figures
resNEtr = plot_lcbd_relationship(richness_NE, lcbd_NE, beta_NE)
resSWtr = plot_lcbd_relationship(richness_SW, lcbd_SW, beta_SW; scaling_value=1000)

# Fix axis dimensions
yticks!(resNEtr[1], 40:2:50)
yticks!(resSWtr[1], 30:2:40)
xticks!(resNEtr[1], -80:4:-60)
xticks!(resSWtr[1], -120:4:-100)
if outcome == "bart"
    local cmax = 450
    resNEtr[2][:clims] = (-Inf, cmax)
    resSWtr[2][:clims] = (-Inf, cmax)
elseif outcome == "raw"
    local cmax = 285
    resNEtr[2][:clims] = (-Inf, cmax)
    resSWtr[2][:clims] = (-Inf, cmax)
end

# Combine subarea figures in 4 panel figure
combined_plot = plot(
    resNEtr,
    resSWtr;
    layout=grid(2, 1),
    size=(900, 600),
    bottommargin=0.0mm,
    title=["a) Northeast subregion" "" "b) Southwest subregion" ""],
    dpi=200,
)

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(combined_plot, joinpath("fig", outcome, "05_$(outcome)_subareas.png"))
end

## Repeat for different subareas

function plot_subareas(
    coords,
    initial_distributions;
    display_coords=coords,
    transform=true,
    kw...
)
    # Get analysis values
    ds = [d[coords] for d in initial_distributions]
    y = Ymatrix(ds)
    r = richness(y, ds[1])
    l = lcbd(y, ds[1]; transform=transform)
    bt = beta_total(y)

    # Plot subareas
    p = plot_lcbd_relationship(r, l, bt; kw...)

    return p
end

# Set coordinates for initial subarea
left = -71.0; right = -64.0; bottom = 46.0; top = 50.0;
coords_subarea = (left=left, right=right, bottom=bottom, top=top)

# Test plot initial subarea
p = plot_subareas(
    coords_subarea,
    distributions;
    formatter=[f -> "$(Int(f))" :plain],
    title=["LCBD" "Relationship"]
)

## Create expanding GIF

# Set initial coordinates
left = -71.0; right = -64.0; bottom = 46.0; top = 50.0;
coords_subarea = (left=left, right=right, bottom=bottom, top=top)

# Set values to control spatial expansion
asp_ratio = 92.60 / 60.75
dim_ratio = (top - bottom)/(right - left)

# Get increasing subareas coordinates
subarea_coords = []
@time while left > -145.0 + asp_ratio && bottom > 20.0 + asp_ratio * dim_ratio
    global left -= asp_ratio
    global bottom -= asp_ratio * dim_ratio
    bbox = (left=left, right=right, bottom=bottom, top=top)
    push!(subarea_coords, bbox)
end

# Plot all expanding subareas
subarea_plots = []
for sc in subarea_coords
    local p = plot_subareas(
        sc,
        distributions;
        formatter=[f -> "$(Int(round(f, digits=0)))" :plain],
        title=["LCBD" "Relationship"],
        dpi=200
    )
    push!(subarea_plots, p)
end

# Create GIF
anim = @animate for p in subarea_plots[Not(1)]
    plot(p)
end
gif(anim, fps=3)

# Export GIF
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    gif(anim, joinpath("fig", outcome, "05_$(outcome)_scaling.gif"); fps=3)
end

## 3 extents comparison

# Extract LCBD & relationship subplots for first, middle & last GIF plots
mid_ind = median(1:length(subarea_plots)) |> round |> Int64
ps = subarea_plots[[1, mid_ind, end]]

# Combine 3 scales in single figure
p = plot(
    deepcopy(ps)...;
    dpi=200,
    layout=(3, 1),
    size=(900, 960),
    title=["a) Regional extent" "" "b) Intermediate extent" "" "c) Continental extent" ""],
    titleloc=:left,
    bottommargin=-2.0mm,
)

# Tweak axes
yticks!(p[3], 34:3:50)
if outcome == "bart"
    yticks!(p[6], 1.5:0.5:5.0)
elseif outcome == "raw"
    yticks!(p[4], 0.5:0.5:2.5)
    yticks!(p[6], 2:1:9)
end
p

# Export figure
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(p, joinpath("fig", outcome, "05_$(outcome)_extents.png"))
end


## Scaling medians figure

# Create empty elements
richness_medians = []
lcbd_medians = []
lcbd_mins = []
lcbd_maxs = []
lcbd_abs_medians = []
beta_values = []
gamma_values = []

# Get analysis values for all subareas
for sc in subarea_coords
    local ds = [d[sc] for d in distributions]
    local y = Ymatrix(ds)
    local r = richness(y, ds[1])
    local l = lcbd(y, ds[1]; relative=false)
    local l_abs = lcbd(y, ds[1]; relative=false)
    local bt = beta_total(y)
    local g = gamma(y)

    push!(richness_medians, median(r))
    push!(lcbd_medians, median(l))
    push!(lcbd_mins, minimum(l))
    push!(lcbd_maxs, maximum(l))
    push!(lcbd_abs_medians, l_abs)
    push!(beta_values, bt)
    push!(gamma_values, g)
end

# Check values
richness_medians
lcbd_medians
beta_values
gamma_values

# Transform to relative values (to plot all variables on same axis)
medians_df = DataFrame(
    idx=eachindex(richness_medians),
    richness=richness_medians ./ maximum(richness_medians),
    lcbd=lcbd_medians ./ maximum(lcbd_medians),
    beta=beta_values ./ maximum(beta_values),
    gamma=gamma_values ./ maximum(gamma_values)
)

# Plot values across scales (step-by-step plots)
medians_p1 = plot(
    medians_df.idx,
    medians_df.richness;
    label="Median richness",
    lw=2,
    xlabel="Subarea extent",
    ylabel="Subarea value (relative to maximum)",
    legend=:bottomright,
    xticks=:none,
    ylim=(0, 1),
    top_margin=mm
)
medians_p2 = plot!(deepcopy(medians_p1), medians_df.lcbd, label="Median LCBD", lw=2)
medians_p3 = plot!(deepcopy(medians_p2), medians_df.beta, label="Total beta diversity", lw=2)
medians_p4 = plot!(deepcopy(medians_p3), medians_df.gamma, label="Gamma diversity", lw=2)

# Export figure
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    # Final plot
    savefig(medians_p4, joinpath("fig", outcome, "05_$(outcome)_medians.png"))
end
