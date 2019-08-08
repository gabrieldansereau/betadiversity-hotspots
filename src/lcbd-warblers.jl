using Distributed
using JLD2
@time include("required.jl")

## Load predictions for all species
@load "../data/predictions-can.jld2" predictions

## Create matrix Y (site-by-species community data table)
begin
    # Get dimensions
    nsites = prod(size(predictions[1]))
    nspecies = length(predictions)
    # Create Y
    Y = zeros(Int64, (nsites, nspecies))
    # Fill Y with community predictions
    @progress for gc in eachindex(predictions[1].grid) # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], predictions)
        # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
        global Y[gc,:] = .!isnan.(R)
    end
end

## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
# Compute BD statistics
resBD = BD(Y)
# Extract LCBD values
LCBDi = resBD.LCBDi

## Arrange LCBD values as grid
# Create empty grid
t_lcbd = zeros(Float64, size(predictions[1]))
# Scale LCBDi values to maximum value
LCBDi = LCBDi./maximum(LCBDi)
# Fill-in grid
for i in eachindex(t_lcbd)
    # Add LCBD values only if prediction != NaN
    t_lcbd[i] = any(Y[i,:] .> 0) ? LCBDi[i] : NaN # ?: is ternary operator, ~if-else
end
# Create SDMLayer with LCBD values
LCBD = SDMLayer(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
sdm_plot = plotSDM(LCBD, type="lcbd")

## Save result
savefig(sdm_plot, "fig/warblers/lcbd-map-can.pdf")
