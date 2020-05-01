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

## Dummy example

# Dummy set
mat = Array(reshape(1.:12., (3,4)))
pos = Array(reshape(101.:112., (3,4)))
# Window size
wsize = (2,2)
# Number of iterations
rowit = size(mat,1) - wsize[1] + 1
colit = size(mat,2) - wsize[2] + 1
totit = rowit * colit
# Cartesian indices
cartinds = findall(!isnan, mat)
mat[cartinds[1]]
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
newmats = get_windows(mat, pos, wsize)

# Reduce to minimum value
@time reduce(min, map(vec, newmats)) |> x-> reshape(x, size(mat))
@time mapreduce(vec, min, newmats) |> x-> reshape(x, size(mat))
# Combine values
@time mapreduce(vec, hcat, newmats)
@time reduce(hcat, map(vec, newmats))
hcat(tmp...)

# Custom mean function
function mean_nonan(x, y)
  mat_combined = mapreduce(vec, hcat, [x, y])
  mat_mean = map(x -> mean(filter(!isnan, x)), eachrow(mat_combined))
  mat_mean = reshape(mat_mean, size(x))
end
# Reduce by mean
mean_nonan(newmats[1], newmats[2])
reduce(mean_nonan, newmats)

## Real example

coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
lonfull = longitudes(distributions[1])
lonNE = longitudes(distributions[1][coords_NE])
indexin(lonNE, lonfull)

distributions_NE = [d[coords_NE] for d in distributions]
wsize = size(distributions_NE[1])

dpos = [CartesianIndices(d.grid) for d in distributions]

@time dwinds = @showprogress map((d, pos) -> get_windows(d.grid, pos, wsize; step = 40), distributions, dpos);

[d[1] for d in dwinds]
