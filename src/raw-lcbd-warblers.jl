using Distributed
using JLD2
@time include("required.jl")

## Load presence-absence data for all species
@load "../data/pres-abs-ebd.jld2" pres_abs
## Load matrix Y
@load "../data/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

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
