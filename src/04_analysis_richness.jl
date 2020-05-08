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
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

#### Species richness
Y = calculate_Y(distributions)
richness = calculate_richness(Y, distributions[1])

## Plot results
richness_plot = plotSDM(richness, c=:viridis,
                        title = "Richness ($outcome distributions)",
                        clim=(0.0, Inf),
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
