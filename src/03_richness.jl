import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

#=
outcome = "sdm"
=#

## Load distribution data for all species
@load "data/jld2/$(outcome)-distributions.jld2" distributions

## Load matrix Y
@load "data/jld2/$(outcome)-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums[inds_notobs] .= NaN
# Reshape to grid format
sums = reshape(sums, size(distributions[1]))

## Create SimpleSDMLayer
richness = SimpleSDMResponse(sums, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

## Plot results
richness_plot = plotSDM(richness, c=:viridis)
heatmap!(richness_plot,
          title = "Richness ($outcome distributions)",
          clim=(0.0, 60.0),
          colorbar_title = "Number of species per site",
          dpi=300)

## Save result
#=
savefig(richness_plot, "fig/$(outcome)/03_$(outcome)_richness.pdf")
=#