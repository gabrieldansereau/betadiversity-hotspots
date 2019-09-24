using Distributed
using JLD2
@time include("../required.jl")

## Load predictions for all species
@load "data/jld2/sdm-predictions.jld2" predictions

## Load matrix Y
@load "data/jld2/sdm-Y-matrices.jld2" Y Ypred Yprob Ytransf inds_pred inds_notpred

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums[inds_notpred] .= NaN
# Reshape to grid format
sums = reshape(sums, size(predictions[1]))

## Create SimpleSDMLayer
richness = SimpleSDMResponse(sums, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
richness_plot = plotSDM(richness, type="lcbd")
title!(richness_plot, "Number of species per site (SDM predictions)")

## Save result
#=
savefig(richness_plot, "fig/sdm/03_sdm_richness.pdf")
=#
