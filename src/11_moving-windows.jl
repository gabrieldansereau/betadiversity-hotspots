import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
outcome = "rf"
# save_figures = true

# Make sure "outcome" is defined
outcome = "rf"
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw', 'sdm' or 'rf'"
elseif !(outcome in ["raw", "sdm", "rf"])
    @warn "'outcome' invalid, must be either 'raw', 'sdm' or 'rf'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Prepare data
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

Ymats = calculate_Ymatrix(distributions)
richness = calculate_richness(Ymats.Y, Ymats.inds_notobs, distributions)
lcbd = calculate_lcbd(Ymats.Yobs, Ymats.Ytransf, Ymats.inds_obs, distributions)

plot(richness, c = :viridis)

## Prepare subareas

# Define subareas
coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
coords_SW = (left = -120.0, right = -100.0, bottom = 30.0, top = 40.0)
distributions_NE = [d[coords_NE] for d in distributions]
# Define window size
wsize = size(distributions_NE[1])

## Moving windows

function get_windows_indices(index_mat, wsize, steps::Tuple{Int64, Int64})
    lat_step, lon_step = steps
    # Get all windows
    winds = []
    for j in 1:lon_step:size(index_mat,2)-(wsize[2]-1), i in 1:lat_step:size(index_mat,1)-(wsize[1]-1)
        subrows, subcols = i:i+(wsize[1]-1), j:j+(wsize[2]-1)
        subinds = index_mat[subrows, subcols]
        push!(winds, subinds)
    end
    return winds
end
get_windows_indices(index_mat, wsize, steps::Int64) = get_windows_indices(index_mat, wsize, (steps, steps))

# Apply moving windows
@time begin
    # Get matrix Y
    Y = calculate_Y(distributions)
    # Create matrix of indices
    index_mat = reshape(eachindex(distributions[1].grid), size(distributions[1])) |> Array
    # Get windows indices
    steps = Int64.(round.(1.0./stride(distributions[1])))
    windows_inds = get_windows_indices(index_mat, wsize, steps)
    # Extract windows from Y
    Ywindows = [Y[vec(winds),:] for winds in windows_inds]
    # Remove ones without observation
    allnan = map(x -> all(isnan, x), Ywindows)
    deleteat!(windows_inds, allnan)
    deleteat!(Ywindows, allnan)
    # Calculate LCBD values for Ywindows
    LCBDwindows = @showprogress [calculate_lcbd(Y, distributions_NE; relative = true) for Y in Ywindows]

    # Expand LCBD windows to match full extent
    LCBDbatch = []
    for i in eachindex(LCBDwindows)
        mat = fill(NaN, size(distributions[1]))
        mat[windows_inds[i]] = LCBDwindows[i][2].grid
        push!(LCBDbatch, mat)
    end
    # Stack windows in single matrix
    LCBDmat = reduce(hcat, map(vec, LCBDbatch));
    # Get mean LCBD value per site
    LCBDmean = map(x -> mean(filter(!isnan, x)), eachrow(LCBDmat))
    # Arrange mean values as layer
    LCBDwindow = similar(distributions[1])
    LCBDwindow.grid .= NaN
    LCBDwindow.grid = reshape(LCBDmean, size(LCBDwindow.grid))
end;
plot(LCBDwindow, c = :viridis)

## Visualize result
function plot_lcbd_windows(richness, lcbd; title = "", kw...)
    p1 = plot(lcbd, c = :viridis, title = "LCBD", colorbar_title = "Relative LCBD score", clim = (0,1))
    p2 = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
                xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                xlim = (1, 45), ylim = (0.0, 1.0))
    if title != ""
        l = @layout [t{.01h}; grid(1,2)]
        ptitle = plot(annotation = (0.5, 0.5, "$title"), framestyle = :none)
        p = plot(ptitle, p1, p2, layout = l; kw...)
    else
        l = @layout [a b]
        p = plot(p1, p2, layout = l; kw...)
    end
    return p
end
window_full = plot_lcbd_windows(richness, LCBDwindow, dpi = 150, size = (900,300))
window_NE = plot_lcbd_windows(richness[coords_NE], LCBDwindow[coords_NE], dpi = 150, size = (900,300))
window_SW = plot_lcbd_windows(richness[coords_SW], LCBDwindow[coords_SW], dpi = 150, size = (900,300))
ptitle_NE = plot(annotation = (0.5, 0.5, "Northeast subarea", 16), framestyle = :none)
ptitle_SW = plot(annotation = (0.5, 0.5, "Southwest subarea", 16), framestyle = :none)
l = @layout [t1{.01h}; p1;
             t2{.01h}; p2]
psubareas = plot(ptitle_NE, window_NE, ptitle_SW, window_SW, layout = l,
                 size = (900,600), dpi = 150)

# Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(window_full, joinpath("fig", outcome, "11_$(outcome)_moving-windows_full.png"))
    savefig(psubareas, joinpath("fig", outcome, "11_$(outcome)_moving-windows_subareas.png"))
end

## GIF
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
    p = plot_lcbd_windows(richness[coords_subarea], LCBDwindow[coords_subarea],
                            formatter = f -> "$(round(f, digits = 1))",
                            aspect_ratio = [asp_ratio :auto],
                            dpi = 150, size = (900,300))
    push!(subarea_plots, p)
end
anim = @animate for p in subarea_plots[Not(1)]
    plot(p)
end
gif(anim, fps = 3)
if (@isdefined save_figures) && save_figures == true
    gif(anim, joinpath("fig", outcome, "11_$(outcome)_moving-windows.gif"), fps = 3)
end

## 3 scales comparison

# Extract LCBD & relationship subplots for first, middle, last GIF plots
ps = []
mid_ind = median(1:length(subarea_plots)) |> round |> Int64
for p in subarea_plots[[1, mid_ind, end]]
    p_lcbd = p[1][1][:plot_object]
    p_rel = p[2][1][:plot_object]
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
    savefig(p, joinpath("fig/", outcome, "11_$(outcome)_moving-windows_3scales.png"))
end
