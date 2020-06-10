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

## Get full-scale comparisons
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
if (@isdefined subset_qc) && subset_qc == true
    coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
    distributions = [d[coords_qc] for d in distributions]
end
Ysdm = calculate_Y(distributions)
richness_sdm = calculate_richness(Ysdm, distributions[1])
lcbd_sdm = calculate_lcbd(Ysdm, distributions[1])

## Arrange predictions as layers
richness_pred = similar(richness_sdm)
richness_pred.grid[:] = predictions[:, 1]
lcbd_pred = similar(lcbd_sdm)
lcbd_pred.grid[:] = predictions[:, 2]

## Plot predictions
richness_plot = plotSDM(richness_pred, c = :viridis,
                        title = "Direct richness predictions ($(uppercase(outcome)))",
                        colorbar_title = "Predicted number of species",
                        )
lcbd_plot = plotSDM(lcbd_pred, c = :viridis,
                    title = "Direct LCBD predictions ($(uppercase(outcome)))",
                    colorbar_title = "LCBD scores",
                    clim = (0,1),
                    )

## Plot differences

# Custom function
function plot_difference(prediction, comparison; absolute = false, hist = false, kw...)
    # Calculate difference
    difference = similar(prediction)
    difference.grid = prediction.grid .- comparison.grid
    # Get absolute values (optional)
    if absolute
        replace!(x -> !isnan(x) ? abs(x) : x, difference.grid)
    end

    # Plot difference
    if hist # as a histogram (optional)
        diff_hist = histogram(filter(!isnan, difference.grid); kw...)
        return diff_hist
    else # as a heatmap (default)
        diff_plot = plotSDM(difference; kw...)
        return diff_plot
    end
end

# Richness difference
richness_diff = plot_difference(richness_pred, richness_sdm;
                            c = :diverging, clim = (-30, 30),
                            title = "Direct richness predictions difference ($(uppercase(outcome)))",
                            colorbar_title = "Difference from SDM-predicted richness",
                            )
richness_absdiff = plot_difference(richness_pred, richness_sdm;
                               absolute = true,
                               c = :inferno, clim = (-Inf, Inf),
                               title = "Direct richness predictions difference ($(uppercase(outcome)))",
                               colorbar_title = "Difference from SDM-predicted richness (absolute)",
                               )

# LCBD difference
lcbd_diff = plot_difference(lcbd_pred, lcbd_sdm,
                            absolute = true,
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
