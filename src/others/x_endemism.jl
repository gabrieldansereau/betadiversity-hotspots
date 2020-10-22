import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

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

## Load distributions for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions spenames speindex
## Load matrix Y
Y = calculate_Y(distributions)

## Weighted endemism
# AOO (Area of occurrence): Species add 1/(number of sites occupied) to site scores

# Number of occupied sites per species
n_occ = map(x -> sum(filter(!isnothing, x)), eachcol(Y))
# Weight site score by number of occupied sites per species
Yweight = hcat([replace(Y[:,i], 1 => 1/n_occ[i]) for i in 1:size(Y,2)]...)
# Sum scores per site for all species
endemism_scores = fill(nothing, size(Yweight, 1)) |> Array{Union{Nothing, Float32}}
inds_obs = _indsobs(Yweight)
endemism_scores[inds_obs] = sum.(eachrow(Yweight[inds_obs]))
# Check stats
describe(filter(!isnothing, endemism_scores))
sort(filter(!isnothing, endemism_scores), rev=true) # possibly some outliers, scaling problems for visualization
# Reshape to grid format
endemism_grid = reshape(endemism_scores, size(distributions[1]))

# Create SimpleSDMLayer with endemism values
endemism = SimpleSDMResponse(Array(endemism_grid), distributions[1])

## Plot results
# Raw endemism scores
endemism_plot = plotSDM2(endemism, c=:viridis,
                        title = "Endemism ($outcome distributions)",
                        colorbar_title = "Weigthed endemism (area of occurrence)",
                        )
# Quantile endemism scores
endemism_qplot = plotSDM2(quantiles(endemism), c=:viridis,
                         title = "Endemism quantiles ($outcome distributions)",
                         colorbar_title = "Weighted endemism quantile (area of occurrence)",
                         )

## Save result
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) lcbd)"
    savefig(endemism_plot, joinpath("fig", outcome, "07_$(outcome)_endemism.png"))
else
    @info "Figures not saved ($(outcome) lcbd)"
end
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) lcbd)"
    savefig(endemism_qplot, joinpath("fig", "quantiles", "07_$(outcome)_endemism_quantiles.png"))
else
    @info "Figures not saved ($(outcome) lcbd)"
end
