using Distributed
using JLD2
@time include("required.jl")

## Load presence-absence data for all species
@load "../data/pres-abs-ebd.jld2" pres_abs

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
# Export matrix Y
@save "../data/Y-pres-abs-ebd.jld2" Y
@load "../data/Y-pres-abs-ebd.jld2" Y

# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
# Select sites with predictions only
Ypred = Y[inds_pred,:]

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

## Data transformation (using vegan in R)
using RCall
@rput Ypred
begin
    R"""
        library(vegan)
        Ytransf <- decostand(Ypred, "hel")
    """
end
@rget Ytransf

## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
## Option 2: Calculate LCBD only for sites with predictions
# Compute BD statistics
resBDpred = BD(Ypred)
resBDtransf = BD(Ytransf)
# Extract LCBD values
LCBDi = resBDpred.LCBDi
LCBDi_tr = resBDtransf.LCBDi
# Scale LCBDi values to maximum value
LCBDi = LCBDi./maximum(LCBDi)
LCBDi_tr = LCBDi_tr./maximum(LCBDi_tr)

## Arrange LCBD values as grid
# Create empty grid
t_lcbd = fill(NaN, size(pres_abs[1]))
t_lcbd_tr = fill(NaN, size(pres_abs[1]))
# Fill in grid
t_lcbd[inds_pred] = LCBDi
t_lcbd_tr[inds_pred] = LCBDi_tr
# Create SDMLayer with LCBD values
LCBD = SDMLayer(t_lcbd, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)
LCBDtr = SDMLayer(t_lcbd_tr, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
lcbd_plot = plotSDM(LCBD, type="lcbd")
title!(lcbd_plot, "LCBD values per site (relative to maximum, raw data)")
lcbdtr_plot = plotSDM(LCBDtr, type="lcbd")
title!(lcbdtr_plot, "LCBD values per site (relative to maximum, hellinger transformed)")

## Save result
#=
savefig(lcbd_plot, "fig/raw-lcbd.pdf")
savefig(lcbdtr_plot, "fig/raw-lcbd-transf.pdf")
=#
