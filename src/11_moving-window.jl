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

## Prepare data
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

distributions[1]

Ymats = calculate_Ymatrix(distributions)
richness = calculate_richness(Ymats.Y, Ymats.inds_notobs, distributions)
lcbd = calculate_lcbd(Ymats.Yobs, Ymats.Ytransf, Ymats.inds_obs, distributions)

plotSDM(richness, c = :viridis)

## Alternative analysis functions

# Previous workflow
@time begin
  Ymats = calculate_Ymatrix(distributions);
  @time richness = calculate_richness(Ymats.Y, Ymats.inds_notobs, distributions);
  @time lcbd = calculate_lcbd(Ymats.Yobs, Ymats.Ytransf, Ymats.inds_obs, distributions);
end;

# New workflow
@time begin
  Y = calculate_Y(distributions, transform = false);
  @time richness = calculate_richness(Y, distributions);
  @time lcbd = calculate_lcbd(Y, distributions; relative = true);
end;

## Moving windows

function get_windows(mat, pos, wsize; step = 1)
  # Get all windows
  newmats = []
  for j in 1:step:size(mat,2)-(wsize[2]-1), i in 1:step:size(mat,1)-(wsize[1]-1)
    subrows, subcols = i:i+(wsize[1]-1), j:j+(wsize[2]-1)
    submat = mat[subrows, subcols]
    subpos = pos[subrows, subcols]
    subinds = indexin(vec(subpos), vec(pos))
    newmat = fill(NaN, size(mat))
    newmat[subinds] = submat
    push!(newmats, newmat)
  end
  return newmats
end

# Define subareas
coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
coords_SW = (left = -120.0, right = -100.0, bottom = 30.0, top = 40.0)
distributions_NE = [d[coords_NE] for d in distributions]
# Smaller test area
distributions_full = copy(distributions)
distributions = copy(distributions_NE)
coords_mini = (left = -80.0, right = -70.0, bottom = 40.0, top = 45.0)
distributions_mini = [d[coords_mini] for d in distributions]
wsize = size(distributions_mini[1])
grid_pos = CartesianIndices(distributions[1].grid)
# Define window size
wsize = size(distributions_NE[1])
# Extract grid positions
grid_pos = CartesianIndices(distributions[1].grid)

## New version
function get_windows_indices(index_mat, wsize; step = 1)
  # Get all windows
  winds = []
  for j in 1:step:size(index_mat,2)-(wsize[2]-1), i in 1:step:size(index_mat,1)-(wsize[1]-1)
    subrows, subcols = i:i+(wsize[1]-1), j:j+(wsize[2]-1)
    subinds = index_mat[subrows, subcols]
    push!(winds, subinds)
  end
  return winds
end
@time begin
  # Get matrix Y
  Y = calculate_Y(distributions)
  # Create matrix of indices
  index_mat = reshape(eachindex(distributions[1].grid), size(distributions[1])) |> Array
  # Get windows indices
  windows_inds = get_windows_indices(index_mat, wsize; step = 10)
  # Extract windows from Y
  Ywindows = [Y[vec(winds),:] for winds in windows_inds]
  # Calculate LCBD values for Ywindows
  LCBDwindows = [calculate_lcbd(Y, distributions_mini; relative = true) for Y in Ywindows]

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
  LCBDwindow.grid = reshape(LCBDmean, size(LCBDwindow.grid))
  LCBDwindow
end;
LCBDwindow_new = copy(LCBDwindow)
plotSDM(LCBDwindow_new, c = :viridis)

## Previous option
@time begin
  # Get distribution windows (per species)
  dwindows = map(d -> get_windows(d.grid, grid_pos, wsize; step = 10), distributions);

  # Arrange windows in a 2D array
  dmats = [d[i] for i in eachindex(dwindows[1]), d in dwindows]
  size(dmats)
  # Arrange as layers in 2D array
  dlayers = [SimpleSDMResponse(d, distributions[1].left, distributions[1].right,
                                distributions[1].bottom, distributions[1].top) for d in dmats]
  size(dlayers)

  # Get Ymatrices
  Ybatch = [calculate_Y(distributions) for distributions in eachrow(dlayers)]
  # Remove ones without observation
  filter!(x -> !all(isnan, x), Ybatch)
  # Get LCBD values
  LCBDbatch = [calculate_lcbd(Y, distributions; relative = true) for Y in Ybatch]

  # Arrange LCBD values as matrix
  LCBDmat = reduce(hcat, [vec(LCBD[2].grid) for LCBD in LCBDbatch])
  # Get mean LCBD value per site
  LCBDmean = map(x -> mean(filter(!isnan, x)), eachrow(LCBDmat))
  # Arrange mean values as layer
  LCBDwindow = similar(distributions[1])
  LCBDwindow.grid = reshape(LCBDmean, size(LCBDwindow.grid))
end;
LCBDwindow_prev = copy(LCBDwindow)
plotSDM(LCBDwindow_prev, c = :viridis)

# Time comparison
# prev: 1.7 sec, 5.50M alloc, 1.662 GB
# new: 0.5sec, 3.07M alloc, 334.514 MB
# Compare number of elements in options
nprev = (62*7000) + (62*28*7000) + (28*62*7000) + (28*7000) + (7000)
nnew = (62*7000) + (1*62*7000) + (28*62*1800) + (28*1800) + (28*7000) + (7000)
nnew/nprev # 17% of previous nvals
# Compare results
filter(!isnan, LCBDwindow_prev.grid) == filter(!isnan, LCBDwindow_new.grid)

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

# Export figures
savefig(window_full, joinpath("fig", outcome, "11_$(outcome)_moving-window_full.png"))
savefig(window_NE, joinpath("fig", outcome, "11_$(outcome)_moving-window_NE.png"))
savefig(window_SW, joinpath("fig", outcome, "11_$(outcome)_moving-window_SW.png"))
