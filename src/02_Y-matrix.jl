import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "sdm" # desired outcome (required)
# save_data = true # should data files be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw' or 'sdm'"
elseif (outcome != "raw" && outcome != "sdm")
    @warn "'outcome' invalid, must be either 'raw' or 'sdm'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Load presence-absence data for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions spenames speindex

## Create matrix Y (site-by-species community data table)
Y = calculate_Y(distributions)

## Export results
# save_data = true # should data files be overwritten (optional)
if (@isdefined save_data) && save_data == true
    # Export data
    @info "Data exported to file ($(outcome) Y matrix)"
    @save joinpath("data", "jld2", "$(outcome)-Y-matrices.jld2") Y
else
    # Load data
    @info "Data imported from file ($(outcome) Y matrix)"
    @load joinpath("data", "jld2", "$(outcome)-Y-matrices.jld2") Y
end

## Visualize patterns in Y
# Heatmap of Y
Yobs = _Yobs(Y)
heat_y = heatmap(Yobs, title = "$(titlecase(outcome)) matrix Y (unsorted)",
                   ylabel = "Site number", xlabel = "Species number")
# Sort Y by rows & columns sums
rowsum = sum.(eachrow(Yobs))
colsum = sum.(eachcol(Yobs))
sortedrows = sortperm(rowsum)
sortedcols = sortperm(colsum, rev=true)
heat_sortrowcol = heatmap(Yobs[sortedrows, sortedcols],
                          title = "$(titlecase(outcome)) matrix Y (sorted by row and column sums)",
                          ylabel = "Site number", xlabel = "Species number", dpi=300)

## Export results
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) Y matrix)"
    savefig(heat_sortrowcol, joinpath("fig", outcome, "02_$(outcome)_Y-rowcolsorted.png"))
else
    @info "Figures not saved ($(outcome) Y matrix)"
end

#= # Funny looking smudge 😛
heatmap(sort(Yobs, dims=1, by=sum))
=#
