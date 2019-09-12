using Distributed
using JLD2
@time include("required.jl")

## Load predictions for all species
@load "../data/predictions-am-larger2.jld2" predictions

## Load matrix Y
@load "../data/sdm-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

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
t_lcbd = fill(NaN, size(predictions[1]))
# Fill in grid
t_lcbd[inds_pred] = LCBDi
# Create SDMLayer with LCBD values
LCBD = SDMLayer(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot results
lcbd_plot = plotSDM(LCBD, type="lcbd")
title!(lcbd_plot, "LCBD values per site (relative to maximum)")

## Save result
# savefig(lcbd_plot, "fig/sdm-lcbd.pdf")
