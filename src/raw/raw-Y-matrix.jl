using Distributed
using JLD2
@time include("required.jl")

## Load presence-absence data for all species
@load "data/jld2/pres-abs-ebd.jld2" pres_abs

## Create matrix Y (site-by-species community data table)
begin
    # Get dimensions
    nsites = prod(size(pres_abs[1]))
    nspecies = length(pres_abs)
    # Create Y
    Y = zeros(Int64, (nsites, nspecies))
    # Fill Y with community predictions
    @progress for gc in 1:nsites # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], pres_abs)
        # Fill Y with binary values
        global Y[gc,:] = isone.(R)
    end
end

## Create matrix Ypred (only sites with observations)
# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
inds_notpred = findall(.!sites_pred)
# Select sites with predictions only
Ypred = Y[inds_pred,:]

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
@save "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred
# Test import
@load "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

## Visualize patterns in Y
# Heatmap of Y
heat_raw = heatmap(Ypred, title = "Presence-absence matrix Y (unsorted)",
                   ylabel = "Site number", xlabel = "Species number")
# Sort Y by rowsum (I think?)
Ysort = sortslices(Ypred, dims = 1)
# Heatmap of Y_sort
heat_sort = heatmap(Ysort, title = "Presence-absence matrix Y (sortslices)",
                   ylabel = "Site number", xlabel = "Species number")
# Check order
sum(eachcol(Ysort))
sum(eachcol(Ysort))[end-10:end] # not quite in order...
# Custom sorting
rowsum = sum.(eachrow(Ypred))
colsum = sum.(eachcol(Ypred))
sortedrows = sortperm(rowsum)
sortedcols = sortperm(colsum, rev=true)
heat_sortrow = heatmap(Ypred[sortedrows, :], title = "Presence-absence matrix Y (sorted by row sums)",
                   ylabel = "Site number", xlabel = "Species number")
heat_sortrowcol = heatmap(Ypred[sortedrows, sortedcols], title = "Presence-absence matrix Y (sorted by row and column sums)",
                   ylabel = "Site number", xlabel = "Species number")
# Export results
#=
savefig(heat_raw, "fig/raw-Y-unsorted.png")
savefig(heat_sortrow, "fig/raw-Y-rowsorted.png")
savefig(heat_sortrowcol, "fig/raw-Y-rowcolsorted.png")
=#
#= # Funny looking smudge 😛
heatmap(sort(Ypred, dims=1, by=sum))
=#
