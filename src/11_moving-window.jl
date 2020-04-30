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

distributions[1]

Ymats = calculate_Ymatrix(distributions)
richness = calculate_richness(Ymats.Y, Ymats.inds_notobs, distributions)
lcbd = calculate_lcbd(Ymats.Yobs, Ymats.Ytransf, Ymats.inds_obs, distributions)

plotSDM(richness, c = :viridis)

eachindex(richness.grid)

mat = Array(reshape(1.:12., (3,4)))
pos = Array(reshape(101.:112., (3,4)))

rownum = 2
colnum = 2

rowit = size(mat,1) - rownum + 1
colit = size(mat,2) - colnum + 1
totit = rowit * colit

cartinds = findall(!isnan, mat)
tmp[cartinds[1]]

newmats = []
@time for j in 1:size(mat,2)-(colnum-1), i in 1:size(mat,1)-(rownum-1)
  subrows, subcols = i:i+(rownum-1), j:j+(colnum-1)
  submat = mat[subrows, subcols]
  subpos = pos[subrows, subcols]
  subinds = indexin(vec(subpos), vec(pos))
  newmat = fill(NaN, size(mat))
  newmat[subinds] = submat
  push!(newmats, newmat)
end
newmats

function sum_zat(m1, m2)
    n1 = copy(p1.grid)
    for i in eachindex(p1.grid)
        n1[i] = min(p1.grid[i], p2.grid[i])
    end
    return SimpleSDMResponse(n1, p1.left, p1.right, p1.bottom, p1.top)
end

mapreduce(minimum, newmats)
@time reduce(min, map(vec, newmats)) |> x-> reshape(x, size(mat));
@time mapreduce(vec, min, newmats) |> x-> reshape(x, size(mat));

@time mapreduce(vec, hcat, newmats);
@time reduce(hcat, map(vec, newmats))
hcat(tmp...)

coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
lonfull = longitudes(distributions[1])
lonNE = longitudes(distributions[1][coords_NE])
indexin(lonNE, lonfull)
