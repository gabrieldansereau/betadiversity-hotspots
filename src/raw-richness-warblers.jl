using Distributed
using JLD2
@time include("required.jl")

## Load presence-absence data for all species
@load "../data/pres-abs-ebd.jld2" pres_abs

## Load matrix Y
@load "../data/Y-pres-abs-ebd.jld2" Y

## Get index of sites with observations
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
inds_notpred = findall(.!sites_pred)

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums[inds_notpred] .= NaN
# Reshape to grid format
sums = reshape(sums, size(pres_abs[1]))

## Create SDMLayer
richness = SDMLayer(sums, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
richness_plot = plotSDM(richness, type="lcbd")
title!(richness_plot, "Number of species observed per site")

## Save result
# savefig(richness_plot, "fig/richness-ebd-raw.pdf")

#### Calculate LCBDs (same as in presence-absence.jl)
## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
## Option 2: Calculate LCBD only for sites with predictions
# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
# Select sites with predictions only
Ypred = Y[inds_pred,:]
# Compute BD statistics
resBDpred = BD(Ypred)
# Extract LCBD values
LCBDi = resBDpred.LCBDi
# Scale LCBDi values to maximum value
LCBDi = LCBDi./maximum(LCBDi)

## Arrange LCBD values as grid
# Create empty grid
t_lcbd = fill(NaN, size(pres_abs[1]))
# Fill in grid
t_lcbd[inds_pred] = LCBDi
# Create SDMLayer with LCBD values
LCBD = SDMLayer(t_lcbd, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
lcbd_plot = plotSDM(LCBD, type="lcbd")
title!(lcbd_plot, "LCBD values per site (relative to maximum)")

#### Richness-LCBD relationship
## Scatterplot LCBD ~ richness
relation_plot = scatter(vec(richness.grid), vec(LCBD.grid),
        markersize = 1,
        yticks = 0.0:0.20:1.0,
        title = "Relationship between LCBD and species richness",
        xlabel = "Species richness (number of species)", ylabel = "LCBD (relative to maximum)",
        legend = :none, grid=:none)

## Save result
savefig(relation_plot, "fig/relation-ebd-raw.pdf")
