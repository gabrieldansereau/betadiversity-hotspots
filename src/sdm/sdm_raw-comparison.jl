richnessusing Distributed
using JLD2
@time include("../required.jl")

## Load LCBD results on presence-absence datga
@time include("../raw/03_raw_richness.jl")
@time include("../raw/05_raw_lcbd.jl")

richness_raw = richness
LCBDraw = copy(LCBD)
inds_pred_raw = copy(inds_pred)

## Load LCBD results on SDM predictions
@time include("03_sdm_richness.jl")
@time include("05_sdm_lcbd.jl")

richness_sdm = richness
LCBDsdm = copy(LCBD)
inds_pred_sdm = copy(inds_pred)

#### Compare richness
richness = DataFrame(raw = richness_raw.grid[inds_pred_raw], sdm = richness_sdm.grid[inds_pred_raw])
comparison_richness = scatter(richness.raw, richness.sdm,
                              markersize = 1,
                              title = "Relationship between observed and predicted richness values",
                              xlabel = "Observed richness", ylabel = "Predicted richness",
                              legend = :none, grid = :none)

histogram(richness.raw, title = "Distribution of observed richness",
          xlabel = "Observed richness", ylabel = "Count")
histogram(richness.sdm, title = "Distribution of predicted richness",
          xlabel = "Predicted richness", ylabel = "Count")

richness.diff = richness.sdm - richness.raw
mean(filter(!isnan,richness.diff))
std(filter(!isnan, richness.diff))
variation(filter(!isnan, richness.diff))

histogram(richness.diff, title = "Distribution of richness prediction offset",
        xlabel = "Richness prediction offset", ylabel = "Count")
scatter(richness.raw, richness.diff, markersize = 1,
        title = "Prediction offset according to observed value",
        xlabel = "Observed richness", ylabel = "Prediction offset")


#### Compare LCBDs
## 1. Binary predictions (no transformations)
lcbd_raw = LCBDraw[1].grid[inds_pred_raw]
lcbd_sdm = LCBDsdm[1].grid[inds_pred_raw]
## 2. Transformed raw data, binary predictions
#=
lcbd_raw = LCBDraw[2].grid[inds_pred_raw]
lcbd_sdm = LCBDsdm[1].grid[inds_pred_raw]
=#
## 3. Transformed raw data & predictions
#=
lcbd_raw = LCBDraw[2].grid[inds_pred_raw]
lcbd_sdm = LCBDsdm[2].grid[inds_pred_raw]
=#

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
histogram(lcbd.diff, title = "Distribution of LCBD prediction offset",
          xlabel = "LCBD prediction offset", ylabel = "Count")
