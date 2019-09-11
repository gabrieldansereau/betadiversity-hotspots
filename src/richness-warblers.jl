using Distributed
using JLD2
@time include("required.jl")

## Load predictions for all species
@load "../data/predictions-ebd.jld2" predictions

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

## Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
inds_notpred = findall(.!sites_pred)

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums[inds_notpred] .= NaN
# Reshape to grid format
sums = reshape(sums, size(predictions[1]))

## Calculate zscores for each site
# Mean number if species
meanspecies = mean(filter(!isnan, sums))
# Standard deviation
sdspecies = std(filter(!isnan, sums))
# Absolute z-scores
zscores = abs.(sums .- meanspecies)/sdspecies

## Create SDMLayer
richness = SDMLayer(sums, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)
richness_zscores = SDMLayer(zscores, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
richness_plot = plotSDM(richness, type="lcbd")
title!(richness_plot, "Number of species per site (raw)")
richness_zscores_plot = plotSDM(richness_zscores, type="lcbd")
title!(richness_zscores_plot, "Number of species per site (absolute Z-scores)")

## Save result
# savefig(richness_plot, "fig/sdm-richness.pdf")
# savefig(richness_zscores_plot, "fig/sdm-richness-zscores.pdf")
