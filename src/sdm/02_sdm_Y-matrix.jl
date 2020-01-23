import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

outcome = "sdm"

## Load presence-absence data for all species
@load "data/jld2/$(outcome)-distributions-landcover.jld2" distributions

## Create matrix Y (site-by-species community data table)
begin
    # Get dimensions
    nsites = prod(size(distributions[1]))
    nspecies = length(distributions)
    # Create Y
    Y = zeros(Int64, (nsites, nspecies))
    Yprob = zeros(Float64, (nsites, nspecies))
    # Fill Y with community distributions
    @progress for gc in 1:nsites # loop for all sites
        # Group distributions for all species in site
        R = map(x -> x.grid[gc], distributions)
        # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
        global Y[gc,:] = .!isnan.(R)
        # Fill Yprob with predicted probabilities
        replace!(R, NaN => 0.0)
        global Yprob[gc,:] = R
    end
end

## Create matrix Ypred (only sites with observations)
# Get index of sites with distributions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
inds_notpred = findall(.!sites_pred)
# Select sites with distributions only
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
@save "data/jld2/$(outcome)-Y-matrices-landcover.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred
# Test import
@load "data/jld2/$(outcome)-Y-matrices-landcover.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred

## Visualize patterns in Y
# Heatmap of Y
heat_$(outcome) = heatmap(Ypred, title = "SDM distributions matrix Y (unsorted)",
                   ylabel = "Site number", xlabel = "Species number")
# Custom sorting
rowsum = sum.(eachrow(Ypred))
colsum = sum.(eachcol(Ypred))
sortedrows = sortperm(rowsum)
sortedcols = sortperm(colsum, rev=true)
heat_sortrow = heatmap(Ypred[sortedrows, :], title = "SDM distributions matrix Y (sorted by row sums)",
                   ylabel = "Site number", xlabel = "Species number")
heat_sortrowcol = heatmap(Ypred[sortedrows, sortedcols], title = "SDM distributions matrix Y (sorted by row and column sums)",
                   ylabel = "Site number", xlabel = "Species number")
heat_prob = heatmap(Yprob[sortedrows, sortedcols], title = "SDM distributions matrix Y (probabilities, sorted by row and column sums)",
                   ylabel = "Site number", xlabel = "Species number")

# Export results
#=
savefig(heat_$(outcome), "fig/$(outcome)/02_$(outcome)_Y-unsorted.png")
savefig(heat_sortrow, "fig/$(outcome)/02_$(outcome)_Y-rowsorted.png")
savefig(heat_sortrowcol, "fig/$(outcome)/02_$(outcome)_Y-rowcolsorted.png")
=#
#= # Funny looking smudge 😛
heatmap(sort(Ypred, dims=1, by=sum))
=#
