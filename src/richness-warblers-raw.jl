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
