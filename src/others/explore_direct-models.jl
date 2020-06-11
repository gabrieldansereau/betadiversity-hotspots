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
                        colorbar_title = "Predicted number of species")

## Plot differences

# Get comparisons richness & lcbd
Yraw = calculate_Y(raw_distributions)
Ysdm = calculate_Y(sdm_distributions)
raw, sdm = [(richness = calculate_richness(Y, raw_distributions[1]),
             lcbd = calculate_lcbd(Y, raw_distributions[1])) for Y in (Yraw, Ysdm)]

# Custom calculation function
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

# Custom plot function
function difference_plot(layer, lim; title = "")
    limrange = (-lim, lim)
    diff_map = plotSDM2(layer,
                        c = :diverging, clim = limrange,
                        title = "Richness difference",
                        colorbar_title = "Difference")
    diff_hist = histogram([filter(x -> x > 0, layer.grid), filter(x -> x <= 0, layer.grid)],
                          bins = :rice, c = [:diverging_r :diverging], legend = :none,
                          ylim = maxrange, # xlabel = "Difference",
                          title = "Distribution of difference values", 
                          orientation = :horizontal)
    diff_title = plot(annotation = (0.5, 0.5, "$(title)"), framestyle = :none)
    l = @layout [t{0.01h}; a{0.6w} b{0.38w}]
    diff_plot = plot(diff_title, diff_map, diff_hist, 
                     size = (800, 400), layout = l)
    return diff_plot
end

# Plot differences
plot_pred_sdm = difference_plot(diff_pred_sdm, maxlim, 
                                title = "Direct richness $(uppercase(outcome)) vs SDM richness")
plot_pred_raw = difference_plot(diff_pred_raw, maxlim, 
                                title = "Direct richness $(uppercase(outcome)) vs raw richness")
plot_sdm_raw  = difference_plot(diff_sdm_raw, maxlim, 
                                title = "SDM richness vs raw richness")

## Extras

# Absolute richness difference
richness_absdiff = plotSDM2(difference(pred.richness, sdm.richness; absolute = true),
                            c = :inferno, # clim = (-Inf, Inf),
                            title = "Direct richness predictions difference ($(uppercase(outcome)))",
                            colorbar_title = "Difference from SDM-predicted richness (absolute)")

# LCBD predictions
lcbd_plot = plotSDM2(lcbd_pred, c = :viridis,
                    title = "Direct LCBD predictions ($(uppercase(outcome)))",
                    colorbar_title = "LCBD scores",
                    clim = (0,1),
                    )

# LCBD difference
lcbd_diff = plotSDM2(difference(pred.lcbd, sdm.lcbd; absolute = true),
                     c = :inferno, # clim = (-Inf, Inf),
                     title = "Direct richness predictions difference ($(uppercase(outcome)))",
                     colorbar_title = "Difference in predicted LCBD (absolute)",
                     )

## Export figures
# save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(richness_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_direct-richness.png"))
    savefig(plot(plot_pred_sdm, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-richness_pred-sdm.png"))
    savefig(plot(plot_pred_raw, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-richness_pred-raw.png"))
    savefig(plot(plot_sdm_raw,  dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-richness_sdm-raw.png"))    
end
# save_extras = true
if (@isdefined save_extras) && save_extras == true
    savefig(plot(richness_absdiff_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-richness-abs.png"))
    savefig(plot(lcbd_plot, dpi = 150),     joinpath("fig", outcome, "x_$(outcome)_direct-lcbd.png"))
    savefig(plot(lcbd_diff_plot, dpi = 150), joinpath("fig", outcome, "x_$(outcome)_diff-lcbd.png"))
end