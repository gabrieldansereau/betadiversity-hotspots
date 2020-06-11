import Pkg
Pkg.activate(".")
using RCall
R"source(file.path('src', 'required.R'))"; # bug with `velox` if not called here
using Distributed
@time include(joinpath("..", "required.jl"))

## Conditional arguments
outcome = "bart"
# outcome = "rf"
# save_figures = true
# subset_qc = true

# Subset to QC data (optional)
if (@isdefined subset_qc) && subset_qc == true
    @rput subset_qc
end

## Train & apply models (RF or BART)
model_script = joinpath("src", "others", "explore_direct-models_$(outcome).R")
@rput model_script
@time begin
    R"""
    source(model_script)
    """
end
@rget results

## Fix missing values
predictions = replace(Array(results[:predictions]), missing => NaN)

## Get comparison layers
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
raw_distributions = distributions
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
sdm_distributions = distributions
if (@isdefined subset_qc) && subset_qc == true
    coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
    raw_distributions = [d[coords_qc] for d in raw_distributions]
    sdm_distributions = [d[coords_qc] for d in sdm_distributions]
end

## Arrange predictions as layers
richness_pred = copy(raw_distributions[1])
richness_pred.grid[:] = predictions[:, 1]
lcbd_pred = copy(raw_distributions[1])
lcbd_pred.grid[:] = predictions[:, 2]
pred = (richness = richness_pred, lcbd = lcbd_pred)

## Plot predictions
richness_plot = plotSDM2(richness_pred, c = :viridis,
                        title = "Direct richness predictions ($(uppercase(outcome)))",
                        colorbar_title = "Predicted number of species",
                        )
lcbd_plot = plotSDM2(lcbd_pred, c = :viridis,
                    title = "Direct LCBD predictions ($(uppercase(outcome)))",
                    colorbar_title = "LCBD scores",
                    clim = (0,1),
                    )

## Plot differences

# Get comparisons richness & lcbd
Yraw = calculate_Y(raw_distributions)
Ysdm = calculate_Y(sdm_distributions)
raw, sdm = [(richness = calculate_richness(Y, raw_distributions[1]),
             lcbd = calculate_lcbd(Y, raw_distributions[1])) for Y in (Yraw, Ysdm)]

# Custom function
function difference(prediction, comparison; absolute = false)
    # Calculate difference
    diff_layer = similar(prediction)
    diff_layer.grid = prediction.grid .- comparison.grid
    # Get absolute values (optional)
    if absolute
        replace!(x -> !isnan(x) ? abs(x) : x, diff_layer.grid)
    end
    return diff_layer
end

# Calculate differences
diff_pred_raw = difference(pred.richness, raw.richness)
diff_pred_sdm = difference(pred.richness, sdm.richness)
diff_sdm_raw  = difference(sdm.richness,  raw.richness)
diffs = [diff_pred_raw, diff_pred_sdm, diff_sdm_raw]
map(x -> extrema(filter(!isnan, x.grid)), diffs)
maxlim = 30
maxrange = (-maxlim, maxlim)

# Richness difference (pred-raw)
map_pred_raw = plotSDM2(diff_pred_raw,
                        c = :diverging, clim = maxrange,
                        title = "Richness difference",
                        colorbar_title = "Difference",
                        )
hist_pred_raw = histogram(filter(!isnan, diff_pred_raw.grid), 
                          xlim = maxrange, xlabel = "Difference", legend = :none,
                          title = "Distribution of difference values", bins = 50)
title_pred_raw = plot(annotation = (0.5, 0.5, "Direct richness $(uppercase(outcome)) vs raw richness"), framestyle = :none)
plot_pred_raw = plot(title_pred_raw, map_pred_raw, hist_pred_raw, 
                     size = (800, 400), layout = @layout [t{0.01h}; a{0.6w} b{0.38w}])

# Richness difference (pred-sdm)
map_pred_sdm = plotSDM2(diff_pred_sdm,
                        c = :diverging, clim = maxrange,
                        title = "Richness difference",
                        colorbar_title = "Difference",
                        )
hist_pred_sdm = histogram(filter(!isnan, diff_pred_sdm.grid),
                          xlim = maxrange, xlabel = "Difference", legend = :none,
                          title = "Distribution of difference values", bins = 50)
title_pred_sdm = plot(annotation = (0.5, 0.5, "Direct richness $(uppercase(outcome)) vs SDM richness"), framestyle = :none)
plot_pred_sdm = plot(title_pred_sdm, map_pred_sdm, hist_pred_sdm, 
                     size = (800, 400), layout = @layout [t{0.01h}; a{0.6w} b{0.38w}])

# Richness difference (sdm-raw)
richness_diff = plotSDM2(diff_sdm_raw,
                        c = :diverging, clim = (-30, 30),
                        title = "Richness difference",
                        colorbar_title = "Difference",
                        )
richness_hist = histogram(filter(!isnan, diff_sdm_raw.grid),
                          xlim = maxrange, xlabel = "Difference", legend = :none,
                          title = "Distribution of difference values", bins = 50)
title_sdm_raw = plot(annotation = (0.5, 0.5, "SDM richness vs raw richness"), framestyle = :none)
plot_sdm_raw = plot(title_sdm_raw, map_sdm_raw, hist_sdm_raw, 
                     size = (800, 400), layout = @layout [t{0.01h}; a{0.6w} b{0.38w}])

## Extras

# Absolute richness difference
richness_absdiff = plotSDM(difference(pred.richness, sdm.richness; absolute = true),
                           c = :inferno, clim = (-Inf, Inf),
                           title = "Direct richness predictions difference ($(uppercase(outcome)))",
                           colorbar_title = "Difference from SDM-predicted richness (absolute)",
                           )

# LCBD difference
lcbd_diff = plotSDM(difference(pred.lcbd, sdm.lcbd; absolute = true),
                    c = :inferno, clim = (-Inf, Inf),
                    title = "Direct richness predictions difference ($(uppercase(outcome)))",
                    colorbar_title = "Difference in predicted LCBD (absolute)",
                    )

## Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(richness_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_direct-richness.png"))
    savefig(plot(lcbd_plot, dpi = 150),     joinpath("fig", outcome, "x_$(outcome)_direct-lcbd.png"))

    savefig(plot(richness_diff_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-richness.png"))
    savefig(plot(richness_absdiff_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-richness-abs.png"))
    savefig(plot(lcbd_diff_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-lcbd.png"))
end
