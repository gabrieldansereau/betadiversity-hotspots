import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "sdm" # desired outcome (required)
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
  @warn "'outcome' not defined, must be either 'raw' or 'sdm'"
elseif (outcome != "raw" && outcome != "sdm")
  @warn "'outcome' invalid, must be either 'raw' or 'sdm'"
else
  @info "'outcome' currently set to '$(outcome)'"
end

## Load distribution data for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions spenames speindex

## Load matrix Y
@load joinpath("data", "jld2", "$(outcome)-Y-matrices.jld2") Y Yobs Ytransf inds_obs inds_notobs

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
richness_plot = plotSDM(richness, c=:viridis,
                        title = "Richness ($outcome distributions)",
                        clim=(0.0, 60.0),
                        colorbar_title = "Number of species per site",
                        dpi=300)
richness_qplot = plotSDM(quantiles(richness), c=:viridis,
                         title = "Richness quantiles ($outcome distributions)",
                         colorbar_title = "Number of species per site (quantiles)",
                         dpi=300)


## Save result
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) richness)"
    savefig(richness_plot, joinpath("fig", outcome, "03_$(outcome)_richness.png"))
else
    @info "Figures not saved ($(outcome) richness)"
end
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) richness)"
    savefig(richness_qplot, joinpath("fig", "quantiles", "03_$(outcome)_richness_quantiles.png"))
else
    @info "Figures not saved ($(outcome) richness)"
end
