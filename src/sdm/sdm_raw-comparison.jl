using Distributed
using JLD2
@time include("../required.jl")

## Load LCBD results on presence-absence datga
@time include("../raw/05_raw_lcbd.jl")

LCBDraw = copy(LCBD)
inds_pred_raw = copy(inds_pred)

## Load LCBD results on SDM predictions
@time include("05_sdm_lcbd.jl")

LCBDsdm = copy(LCBD)
inds_pred_sdm = copy(inds_pred)

#### Compare LCBDs
## 1. Binary predictions (no transformations)
lcbd_raw = LCBDraw[1].grid[inds_pred_raw]
lcbd_sdm = LCBDsdm[1].grid[inds_pred_raw]
lcbd = DataFrame(raw = lcbd_raw, sdm = lcbd_sdm)
comparison_plot1 = scatter(lcbd.raw, lcbd.sdm,
        markersize = 1,
        yticks = 0.0:0.20:1.0,
        title = "Relationship between observed and predicted LCBD values",
        xlabel = "Observed LCBD", ylabel = "Predicted LCBD",
        legend = :none, grid=:none)

lcbd.diff = lcbd.sdm - lcbd.raw
mean(filter(!isnan,lcbd.diff))
std(filter(!isnan, lcbd.diff))
variation(filter(!isnan, lcbd.diff))

scatter(lcbd.raw, lcbd.diff, markersize = 1,
        title = "Prediction offset according to observed value",
        xlabel = "Observed LCBD", ylabel = "Prediction offset")

histogram(lcbd.raw, title = "Distribution of observed LCBDs",
          xlabel = "Observed LCBD", ylabel = "Count")
histogram(lcbd.sdm, title = "Distribution of predicted LCBDs",
          xlabel = "Predicted LCBD", ylabel = "Count")
histogram(lcbd.diff, title = "Distribution of prediction offset",
          xlabel = "Prediction offset", ylabel = "Count")
