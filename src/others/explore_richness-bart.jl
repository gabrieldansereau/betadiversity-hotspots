import Pkg
Pkg.activate(".")
using RCall
begin
    R"""
    source(file.path("src", "required.R"))
    """
end
using Distributed
@time include(joinpath("..", "required.jl"))

## Conditional arguments
# save_figures = true
subset_qc = true

## BART model
# Subset to QC data (optional)
if (@isdefined subset_qc) && subset_qc == true
    @rput subset_qc
end
# Train & apply model
@time begin
    R"""
    source(here("src", "others", "explore_richness-bart.R"))
    """
end
@rget pred_df lower_df upper_df

## Fix missing values
predictions = replace(Array(pred_df), missing => NaN)
lower = replace(Array(pred_df), missing => NaN)
upper = replace(Array(pred_df), missing => NaN)

## Get full-scale comparisons
@load joinpath("data", "jld2", "bart-distributions.jld2") distributions
if (@isdefined subset_qc) && subset_qc == true
    coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
    distributions = [d[coords_qc] for d in distributions]
end
Ysdm = calculate_Y(distributions)
richness_sdm = calculate_richness(Ysdm, distributions[1])
lcbd_sdm = calculate_lcbd(Ysdm, distributions[1])

## Arrange predictions as layers
richness_bart = similar(richness_sdm)
richness_bart.grid[:] = predictions[:, 1]
lcbd_bart = similar(lcbd_sdm)
lcbd_bart.grid[:] = predictions[:, 2]

## Plot predictions
richness_plot = plotSDM(richness_bart, c = :viridis,
                        title = "Richness BART predictions",
                        colorbar_title = "Predicted number of species",
                        )
lcbd_plot = plotSDM(lcbd_bart, c = :viridis,
                    title = "LCBD BART predictions",
                    colorbar_title = "LCBD scores",
                    clim = (0,1),
                    )

## Map richness difference
richness_diff = similar(richness_bart)
richness_diff.grid = abs.(richness_bart.grid .- richness_sdm.grid)
richness_diff_plot = plotSDM(richness_diff, c = :inferno, clim = (-Inf, Inf),
                         title = "Predicted richness - BART vs SDM",
                         colorbar_title = "Difference in predicted richness (absolute)",
                         )
histogram(filter(!isnan, richness_diff.grid), bins = 20)

## Map LCBD difference
lcbd_diff = similar(lcbd_bart)
lcbd_diff.grid = abs.(lcbd_bart.grid .- lcbd_sdm.grid)
lcbd_diff_plot = plotSDM(lcbd_diff, c = :inferno, clim = (-Inf, Inf),
                         title = "Predicted LCBD - BART vs SDM",
                         colorbar_title = "Difference in predicted LCBD (absolute)",
                         )
histogram(filter(!isnan, lcbd_diff.grid), bins = 20)

## Export figures
save_figures = true
if (@isdefined save_figures) && save_figures == true
    savefig(plot(richness_plot, dpi = 150), joinpath("fig", "bart", "x_bart_richness-bart.png"))
    savefig(plot(lcbd_plot, dpi = 150),     joinpath("fig", "bart", "x_bart_lcbd-bart.png"))

    savefig(plot(richness_diff_plot, dpi = 150), joinpath("fig", "bart", "x_bart_richness-diff.png"))
    savefig(plot(lcbd_diff_plot, dpi = 150),     joinpath("fig", "bart", "x_bart_lcbd-diff.png"))
end
