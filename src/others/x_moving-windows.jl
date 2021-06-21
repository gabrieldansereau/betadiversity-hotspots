include("required.jl")

## Conditional arguments
# outcome = "rf"
# outcome = "bart"
# save_figures = true

# Make sure "outcome" is defined
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw', 'bio', 'rf', or 'bart'"
elseif !(outcome in ["raw", "bio", "rf", "bart"])
    @warn "'outcome' invalid, must be either 'raw', 'bio', 'rf', or 'bart'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Prepare data
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
spenames = filter(:type => ==("species"), glossary).full_name
distributions = [
    geotiff(
        SimpleSDMPredictor, joinpath("data", "proc", "distributions_$(outcome).tif"), i
    ) for i in eachindex(spenames)
]

Y = calculate_Y(distributions)
richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1]; transform=true)

plot(richness; c=:viridis)

## Prepare subareas

# Define subareas
coords_NE = (left=-80.0, right=-60.0, bottom=40.0, top=50.0)
coords_SW = (left=-120.0, right=-100.0, bottom=30.0, top=40.0)
distributions_NE = [d[coords_NE] for d in distributions]
# Define window size
wsize = size(distributions_NE[1])

## Moving windows

function get_windows_indices(index_mat, wsize, steps::Tuple{Int64,Int64})
    lat_step, lon_step = steps
    # Get all windows
    winds = []
    jrange = 1:lon_step:(size(index_mat, 2) - (wsize[2] - 1))
    irange = 1:lat_step:(size(index_mat, 1) - (wsize[1] - 1))
    for j in jrange, i in irange
        subrows, subcols = i:(i + (wsize[1] - 1)), j:(j + (wsize[2] - 1))
        subinds = index_mat[subrows, subcols]
        push!(winds, subinds)
    end
    return winds
end
function get_windows_indices(index_mat, wsize, steps::Int64)
    return get_windows_indices(index_mat, wsize, (steps, steps))
end

# Apply moving windows
@time begin
    # Get matrix Y
    Y = calculate_Y(distributions)
    # Create matrix of indices
    index_mat = reshape(eachindex(distributions[1].grid), size(distributions[1])) |> Array
    # Get windows indices
    steps = Int64.(round.(1.0 ./ stride(distributions[1])))
    windows_inds = get_windows_indices(index_mat, wsize, steps)
    # Extract windows from Y
    Ywindows = [Y[vec(winds), :] for winds in windows_inds]
    # Remove ones without observation
    allnothing = map(x -> all(isnothing, x), Ywindows)
    deleteat!(windows_inds, allnothing)
    deleteat!(Ywindows, allnothing)
    # Calculate LCBD values for Ywindows
    LCBDwindows = @showprogress [
        calculate_lcbd(Yobs, distributions_NE[1]; transform=true, relative=true) for
        Yobs in Ywindows
    ]

    # Expand LCBD windows to match full extent
    LCBDbatch = []
    for i in eachindex(LCBDwindows)
        mat = fill(nothing, size(distributions[1])) |> Array{Union{Nothing,Float32}}
        mat[windows_inds[i]] = LCBDwindows[i].grid
        push!(LCBDbatch, mat)
    end
    # Stack windows in single matrix
    LCBDmat = reduce(hcat, map(vec, LCBDbatch))
    # Get mean LCBD value per site
    LCBDmean = map(x -> filter(!isnothing, x), eachrow(LCBDmat))
    LCBDmean = map(x -> isempty(x) ? nothing : mean(x), LCBDmean)
    # Arrange mean values as layer
    LCBDwindow = similar(distributions[1])
    LCBDwindow.grid = reshape(LCBDmean, size(LCBDwindow.grid))
end;
plot(LCBDwindow; c=:viridis)

## Visualize result
if outcome == "raw"
    plotfct = :plotSDM2
else
    plotfct = :plot
end
function plot_lcbd_relationship(richness, lcbd; maintitle="", kw...)
    p1 = eval(plotfct)(
        lcbd; c=:viridis, title="LCBD", colorbar_title="Relative LCBD score", clim=(0, 1)
    )
    p2 = histogram2d(
        richness,
        lcbd;
        c=:viridis,
        bins=40,
        title="Relationship",
        xlabel="Richness",
        ylabel="LCBD",
        colorbar_title="Number of sites",
        xlim=(1, 45),
        ylim=(0.0, 1.0),
        bottommargin=3.0mm,
    )
    if maintitle != ""
        l = @layout [t{0.01h}; grid(1, 2)]
        ptitle = plot(; annotation=(0.5, 0.5, "$maintitle"), framestyle=:none)
        p = plot(ptitle, p1, p2; layout=l, size=(900, 300), kw...)
    else
        l = @layout [a b]
        p = plot(p1, p2; layout=l, size=(900, 300), kw...)
    end
    return p
end
window_full = plot_lcbd_relationship(richness, LCBDwindow)
window_NE = plot_lcbd_relationship(
    richness[coords_NE], LCBDwindow[coords_NE]; maintitle="Northeast subarea"
)
window_SW = plot_lcbd_relationship(
    richness[coords_SW], LCBDwindow[coords_SW]; maintitle="Southwest subarea"
)
psubareas = plot(
    window_NE,
    window_SW;
    layout=grid(2, 1),
    size=(900, 600),
    bottommargin=0.0mm,
    title=["" "" "" ""],
)

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(
        plot(window_full; dpi=200),
        joinpath("fig", outcome, "06-0_$(outcome)_moving-windows_full.png"),
    )
    savefig(
        plot(psubareas; dpi=200),
        joinpath("fig", outcome, "06-1_$(outcome)_moving-windows_subareas.png"),
    )
end

## GIF
left = -71.0;
right = -64.0;
bottom = 46.5;
top = 50.0;
dim_ratio = (top - bottom) / (right - left)
asp_ratio = 92.60 / 60.75
coords_subarea = (left=left, right=right, bottom=bottom, top=top)
subarea_plots = []
nplots = 0
@time while left > -145.0 + asp_ratio && bottom > 20.0 + asp_ratio * dim_ratio
    global nplots += 1
    global left -= asp_ratio
    global bottom -= asp_ratio * dim_ratio
    local coords_subarea = (left=left, right=right, bottom=bottom, top=top)
    local p = plot_lcbd_relationship(
        richness[coords_subarea],
        LCBDwindow[coords_subarea];
        formatter=f -> "$(round(f, digits=1))",
        dpi=200,
    )
    push!(subarea_plots, p)
end
anim = @animate for p in subarea_plots[Not(1)]
    plot(p)
end
gif(anim; fps=3)
if (@isdefined save_figures) && save_figures == true
    gif(anim, joinpath("fig", outcome, "06-3_$(outcome)_moving-windows.gif"); fps=3)
end

#### 3 scales comparison

# Extract LCBD & relationship subplots for first, middle, last GIF plots
mid_ind = median(1:length(subarea_plots)) |> round |> Int64
ps = subarea_plots[[1, mid_ind, end]]

# Combine 3 scales
p = plot(
    ps...;
    dpi=200,
    layout=(3, 1),
    size=(900, 900),
    title=["LCBD" "Relationship" "" "" "" ""],
)

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(p, joinpath("fig/", outcome, "06-2_$(outcome)_moving-windows_3scales.png"))
end
