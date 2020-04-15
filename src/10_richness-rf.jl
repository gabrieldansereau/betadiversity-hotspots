import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load distribution data for all species
@load "data/jld2/raw-distributions.jld2" distributions spenames speindex

## Load matrix Y
@load "data/jld2/raw-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

## Richness
function calculate_richness(Y, inds_notobs, distributions)
  ## Get number of species per site
  sums = map(x -> Float64(sum(x)), eachrow(Y))
  # Add NaN for non predicted sites
  sums[inds_notobs] .= NaN
  # Reshape to grid format
  sums = reshape(sums, size(distributions[1]))
  ## Create SimpleSDMLayer
  richness = SimpleSDMResponse(sums, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)
end
richness_NE = calculate_richness(NE.Y, NE.inds_notobs, distributions_NE)
richness_SW = calculate_richness(SW.Y, SW.inds_notobs, distributions_SW)

# Visualize
plot(richness_NE, c = :viridis)
plot(richness_SW, c = :viridis)
