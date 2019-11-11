using Distributed
using JLD2
@time include("../required.jl")

## Load presence-absence data for all species
@load "data/jld2/sdm-predictions.jld2" predictions

## Create matrix Y (site-by-species community data table)
begin
    # Get dimensions
    nsites = prod(size(predictions[1]))
    nspecies = length(predictions)
    # Create Y
    Y = zeros(Int64, (nsites, nspecies))
    Yprob = zeros(Float64, (nsites, nspecies))
    # Fill Y with community predictions
    @progress for gc in 1:nsites # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], predictions)
        # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
        global Y[gc,:] = .!isnan.(R)
        # Fill Yprob with predicted probabilities
        replace!(R, NaN => 0.0)
        global Yprob[gc,:] = R
    end
end

## Create matrix Ypred (only sites with observations)
# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
inds_notpred = findall(.!sites_pred)
# Select sites with predictions only
Ypred = Y[inds_pred,:]
Yprob = Yprob[inds_pred,:]

## Apply Hellinger transformation (using vegan in R)
using RCall
@rput Ypred
begin
    R"""
        library(vegan)
        Ytransf <- decostand(Ypred, "hel")
    """
end
@rget Ytransf

## Export results
# Export matrices & inds_pred (useful to link Y & Ypred-Ytransf)
@save "data/jld2/sdm-Y-matrices.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred
# Test import
@load "data/jld2/sdm-Y-matrices.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred

## Visualize patterns in Y
# Heatmap of Y
heat_sdm = heatmap(Ypred, title = "SDM predictions matrix Y (unsorted)",
                   ylabel = "Site number", xlabel = "Species number")
# Custom sorting
rowsum = sum.(eachrow(Ypred))
colsum = sum.(eachcol(Ypred))
sortedrows = sortperm(rowsum)
sortedcols = sortperm(colsum, rev=true)
heat_sortrow = heatmap(Ypred[sortedrows, :], title = "SDM predictions matrix Y (sorted by row sums)",
                   ylabel = "Site number", xlabel = "Species number")
heat_sortrowcol = heatmap(Ypred[sortedrows, sortedcols], title = "SDM predictions matrix Y (sorted by row and column sums)",
                   ylabel = "Site number", xlabel = "Species number")
heat_prob = heatmap(Yprob[sortedrows, sortedcols], title = "SDM predictions matrix Y (probabilities, sorted by row and column sums)",
                   ylabel = "Site number", xlabel = "Species number")

# Export results
#=
savefig(heat_sdm, "fig/sdm/02_sdm_Y-unsorted.png")
savefig(heat_sortrow, "fig/sdm/02_sdm_Y-rowsorted.png")
savefig(heat_sortrowcol, "fig/sdm/02_sdm_Y-rowcolsorted.png")
=#
#= # Funny looking smudge 😛
heatmap(sort(Ypred, dims=1, by=sum))
=#