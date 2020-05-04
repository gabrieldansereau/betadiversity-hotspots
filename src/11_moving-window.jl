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

# Export figures
savefig(window_full, joinpath("fig", outcome, "11_$(outcome)_moving-window_full.png"))
savefig(window_NE, joinpath("fig", outcome, "11_$(outcome)_moving-window_NE.png"))
savefig(window_SW, joinpath("fig", outcome, "11_$(outcome)_moving-window_SW.png"))
