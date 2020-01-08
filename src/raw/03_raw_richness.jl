Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load presence-absence data for all species
@load "data/jld2/raw-pres-abs.jld2" pres_abs

## Load matrix Y
@load "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums[inds_notpred] .= NaN
# Reshape to grid format
sums = reshape(sums, size(pres_abs[1]))

## Create SimpleSDMLayer
richness = SimpleSDMResponse(sums, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
richness_plot = plotSDM(richness, c=:viridis)
title!(richness_plot, "Number of species observed per site")

## Save result
#=
savefig(richness_plot, "fig/raw/03_raw_richness.pdf")
=#
