using Pkg: Pkg
Pkg.activate(".")
include("required.jl")

## Conditional arguments
# save_additional_figures = true

## Get data
@load joinpath("data", "jld2", "comparison-results.jld2") raw sdm
results = CSV.read(joinpath("data", "proc", "comparison-results.csv"), DataFrame)
residuals_df = CSV.read(joinpath("data", "proc", "comparison-residuals.csv"), DataFrame)

## Plot distributions
# Get colorbar limits
lims_richness = extrema(mapreduce(collect, vcat, [raw.richness, sdm.richness]))
lims_lcbd = extrema(mapreduce(collect, vcat, [raw.lcbd, sdm.lcbd]))

# Plot combined distributions
p1 = plot(raw.richness_plot; clim=lims_richness)
p2 = plot(sdm.richness_plot; clim=lims_richness)
p3 = plotSDM2(
    rescale(raw.lcbd, extrema(raw.lcbd) .* 100_000);
    c=:viridis,
    colorbar_title="LCBD value (x 100,000)",
)
p4 = plotSDM2(
    rescale(sdm.lcbd, extrema(sdm.lcbd) .* 100_000);
    c=:viridis,
    colorbar_title="LCBD value (x 100,000)",
)
combined_plot = plot(
    p1,
    p2,
    p3,
    p4;
    title=[
        "a) Observed richness",
        "b) Predicted richness",
        "c) Observed uniqueness",
        "d) Predicted uniqueness",
    ] |> permutedims,
    titleloc=:left,
    leftmargin=3.0mm,
    size=(940, 680),
    dpi=200,
)

# Save result
if (@isdefined save_additional_figures) && save_additional_figures == true
    savefig(combined_plot, joinpath("fig", "bart", "07_bart_combined-maps.png"))
end

## Difference plots
# Difference between values from SDM models & raw observations
richness_diff = sdm.richness - raw.richness
lcbd_diff = sdm.lcbdcom - raw.lcbdcom

# Custom fonctions for gradients
function SimpleSDMLayers.rescale(x::AbstractArray, m, M)
    return (x .- minimum(x)) ./ (maximum(x) - minimum(x)) .* (M - m) .+ m
end

function rescalegrad(grad, lims; kw...)
    centervalue = abs(lims[1]) / (lims[2] - lims[1])
    rescaled_grad = cgrad(grad, centervalue; kw...)
    return rescaled_grad
end

function subsetgrad(grad, lims; kw...)
    lims_range = range(lims...; length=20)
    subgrad = cgrad([getindex(cgrad(grad; kw...), lims_range)...])
    return subgrad
end

function recentergrad(grad, lims; kw...)
    lmin = minimum(lims)
    lmax = maximum(lims)
    absmax = max(abs(lmin), lmax)
    if abs(lmin) < lmax
        rlims = rescale([-absmax, lmin, absmax], 0, 1)
        subsetgrad(:PuOr, rlims[2:3]; kw...)
    elseif abs(lmin) > lmax
        rlims = rescale([-absmax, lmax, absmax], 0, 1)
        subsetgrad(:PuOr, rlims[1:2]; kw...)
    end
end

# Test functions
lims = extrema(richness_diff)
rescale([lims[1], 0, lims[2]], 0, 1)
rescalegrad(:PuOr, extrema(richness_diff); rev=true)
subsetgrad(:PuOr, (0.2, 1.0); rev=true)
recentergrad(:PuOr, lims; rev=true)

# Custom function to visualize the difference
function difference_plot(layer::T; title="", kw...) where {T<:SimpleSDMLayer}
    # Center colorscale at zero instead of midpoint between extremas
    lims = extrema(layer)
    # centervalue = abs(lims[1])/(lims[2] - lims[1])
    # scalevalues = rescale([lims[1], 0.0, lims[2]], 0, 1)
    diff_map = plotSDM2(
        layer;
        # c = cgrad(:PuOr, centervalue, rev = true), 
        # c = rescalegrad(:PuOr, extrema(layer); rev = true),
        # c = cgrad(:PuOr, scalevalues, rev = true), 
        c=recentergrad(:PuOr, lims; rev=true),
        clims=lims,
        title="Difference map",
        colorbar_title="Difference from observed value",
    )
    diff_hist = histogram(
        layer;
        bins=50,
        c=:PuOr,
        ylims=lims,
        legend=:none,
        title="Difference distribution",
        xlabel="Frequency",
        orientation=:horizontal,
    )
    # diff_title = plot(annotation = (0.5, 0.5, "$(title)"), framestyle = :none)
    # l = @layout [t{0.01h}; a{0.6w} b{0.38w}]
    l = @layout [a{0.6w} b{0.38w}]
    diff_plot = plot(
        # diff_title, 
        diff_map,
        diff_hist,
        size=(850, 340),
        layout=l,
        bottommargin=3.0mm,
        dpi=200,
        # rightmargin = [0mm 5.0mm 0mm], leftmargin = [0mm 5.0mm 5.0mm];
        rightmargin=[5.0mm 0mm],
        leftmargin=[5.0mm 5.0mm],
        kw...,
    )
    return diff_plot
end
richness_diffplot = difference_plot(
    richness_diff;
    # title = "Predicted richness compared to observed richness",
    # yticks = [:auto :auto (-30:10:40, string.(-30:10:40))],
    yticks=[:auto (-30:10:40, string.(-30:10:40))],
    colorbar_title="Richness difference",
    # ylabel = ["" "Latitude" "Richness difference"]
    ylabel=["Latitude" "Richness difference"],
)
lcbd_diff_resc = rescale(lcbd_diff, extrema(lcbd_diff) .* 100_000)
lcbd_diffplot = difference_plot(
    lcbd_diff_resc;
    # title = "Predicted LCBD compared to observed LCBD",
    colorbar_title="LCBD difference (x 100,000)",
    ylabel=["Latitude" "LCBD difference (x 100,000)"],
    # yticks = [:auto :auto (-4:1:3, string.(-4:1:3))],
    yticks=[:auto (-4:1:3, string.(-4:1:3))],
    # ylim = [:auto extrema(latitudes(lcbd_diff_resc)) extrema(lcbd_diff_resc)],
    ylim=[extrema(latitudes(lcbd_diff_resc)) extrema(lcbd_diff_resc)],
    clim=extrema(lcbd_diff_resc),
)
combined_diffplot = plot(
    deepcopy(richness_diffplot),
    deepcopy(lcbd_diffplot);
    layout=(2, 1),
    title=[
        "a) Difference between richness estimates",
        "",
        "b) Difference between LCBD estimates",
        "",
    ] |> permutedims,
    titleloc=:left,
    dpi=200,
    size=(850, 680),
)

# Save figures
if (@isdefined save_additional_figures) && save_additional_figures == true
    savefig(combined_diffplot, joinpath("fig", "bart", "07_bart_comparison-combined.png"))
end

## Residual visualization
# Arrange data
richres_layer = SimpleSDMResponse(
    residuals_df, :richness, similar(raw.richness); latitude=:latitude, longitude=:longitude
)
richres_qp_layer = SimpleSDMResponse(
    residuals_df,
    :richness_qp,
    similar(raw.richness);
    latitude=:latitude,
    longitude=:longitude,
)
richres_nb_layer = SimpleSDMResponse(
    residuals_df,
    :richness_nb,
    similar(raw.richness);
    latitude=:latitude,
    longitude=:longitude,
)
lcbdres_layer = SimpleSDMResponse(
    residuals_df, :lcbd, similar(raw.lcbd); latitude=:latitude, longitude=:longitude
)
lcbdres_br_layer = SimpleSDMResponse(
    residuals_df, :lcbd_br, similar(raw.lcbd); latitude=:latitude, longitude=:longitude
)

# Plot residuals
function residuals_plot(layer::T; title="", kw...) where {T<:SimpleSDMLayer}
    # Center colorscale at zero instead of midpoint between extremas
    lims = extrema(layer)
    diff_map = plotSDM2(
        layer;
        c=recentergrad(:PuOr, lims; rev=true),
        title="Residuals map",
        colorbar_title="Deviance residuals",
        clims=lims,
    )
    diff_hist = histogram(
        layer;
        bins=50,
        c=:PuOr,
        legend=:none,
        title="Residuals distribution",
        xlabel="Frequency",
        ylabel="Deviance residuals",
        orientation=:horizontal,
        ylims=lims,
    )
    # diff_title = plot(annotation = (0.5, 0.5, "$(title)"), framestyle = :none)
    # l = @layout [t{0.01h}; a{0.6w} b{0.38w}]
    l = @layout [a{0.6w} b{0.38w}]
    diff_plot = plot(
        # diff_title, 
        diff_map,
        diff_hist,
        size=(850, 340),
        layout=l,
        bottommargin=3.0mm,
        dpi=200,
        # rightmargin = [0mm 5.0mm 0mm], leftmargin = [0mm 5.0mm 5.0mm];
        rightmargin=[5.0mm 0mm],
        leftmargin=[5.0mm 5.0mm],
        kw...,
    )
    return diff_plot
end
# richness_resplot residuals_plot(richres_layer; title = "Richness Poisson GLM")
richness_qp_resplot = residuals_plot(richres_qp_layer; title="Richness Quasipoisson GLM")
richness_nb_resplot = residuals_plot(
    richres_nb_layer;
    title="Richness Negative Binomial GLM",
    # yticks = [:auto :auto (-3:1:4, string.(-3:1:4))]
    yticks=[:auto (-3:1:4, string.(-3:1:4))],
)
lcbd_resplot = residuals_plot(lcbdres_layer; title="LCBD Gamma GLM")
lcbd_br_resplot = residuals_plot(lcbdres_br_layer; title="LCBD Beta Regression")

combined_resplot = plot(
    deepcopy(richness_nb_resplot),
    deepcopy(lcbd_br_resplot);
    title=[
        "a) Poisson regression residuals between richness estimates",
        "",
        "b) Beta regression residuals between LCBD estimates",
        "",
    ] |> permutedims,
    titleloc=:left,
    layout=(2, 1),
    size=(850, 680),
)

# Save figures
if (@isdefined save_additional_figures) && save_additional_figures == true
    # savefig(richness_resplot, joinpath("fig", "bart", "07_bart_residuals_richness-poisson.png"))
    savefig(combined_resplot, joinpath("fig", "bart", "07_bart_residuals-combined.png"))
end
