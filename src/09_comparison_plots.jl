include("required.jl")

## Conditional arguments
# save_additional_figures = true

## Get data

# Load raw & sdm data
@load joinpath("data", "jld2", "comparison-results.jld2") raw sdm

# Load GLM results
results = CSV.read(joinpath("data", "proc", "comparison-results.csv"), DataFrame)
residuals_df = CSV.read(joinpath("data", "proc", "comparison-residuals.csv"), DataFrame)

## Plot distributions

# Get colorbar limits
lims_richness = extrema(mapreduce(collect, vcat, [raw.richness, sdm.richness]))
lims_lcbd = extrema(mapreduce(collect, vcat, [raw.lcbd, sdm.lcbd]))

# Plot combined distributions
p1 = plot(raw.richness_plot; clim=lims_richness)
p2 = plot(sdm.richness_plot; clim=lims_richness)
p3 = plot_layer(
    rescale(raw.lcbd, extrema(raw.lcbd) .* 100_000);
    c=:viridis,
    colorbar_title="LCBD value (x 100,000)",
)
p4 = plot_layer(
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
    savefig(combined_plot, joinpath("fig", "bart", "09_bart_combined.png"))
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

    # Difference map subpanel
    diff_map = plot_layer(
        layer;
        c=recentergrad(:PuOr, lims; rev=true),
        clims=lims,
        title="Difference map",
        colorbar_title="Difference from observed value",
    )

    # Difference histogram subpanel
    diff_hist = histogram(
        layer;
        bins=50,
        c=:PuOr,
        ylims=lims,
        legend=:none,
        title="Difference distribution",
        xlabel="Relative frequency",
        orientation=:horizontal,
        normalize=:probability
    )

    # Combine subpanels in single plot
    l = @layout [a{0.6w} b{0.38w}]
    diff_plot = plot(
        diff_map,
        diff_hist;
        size=(850, 340),
        layout=l,
        bottommargin=3.0mm,
        dpi=200,
        rightmargin=[5.0mm 0mm],
        leftmargin=[5.0mm 5.0mm],
        kw...,
    )
    return diff_plot
end

# Richness difference plot
richness_diffplot = difference_plot(
    richness_diff;
    yticks=[:auto (-30:10:40, string.(-30:10:40))],
    colorbar_title="Richness difference",
    ylabel=["Latitude" "Richness difference"],
)

# LCBD difference plot
lcbd_diff_resc = rescale(lcbd_diff, extrema(lcbd_diff) .* 100_000)
lcbd_diffplot = difference_plot(
    lcbd_diff_resc;
    colorbar_title="LCBD difference (x 100,000)",
    ylabel=["Latitude" "LCBD difference (x 100,000)"],
    yticks=[:auto (-4:1:3, string.(-4:1:3))],
    ylim=[extrema(latitudes(lcbd_diff_resc)) extrema(lcbd_diff_resc)],
    clim=extrema(lcbd_diff_resc),
)

# Combine richess & LCBD difference plots
combined_diffplot = plot(
    deepcopy(richness_diffplot),
    deepcopy(lcbd_diffplot);
    layout=(2, 1),
    title=[
        "a) Difference between richness estimates" "" "b) Difference between LCBD estimates" ""
    ],
    titleloc=:left,
    dpi=200,
    size=(850, 680),
)
xlims!(combined_diffplot[2], (-0.005, 0.185))
xlims!(combined_diffplot[4], (-0.005, 0.185))

# Save figures
# save_additional_figures = true
if (@isdefined save_additional_figures) && save_additional_figures == true
    savefig(combined_diffplot, joinpath("fig", "bart", "09_bart_difference.png"))
end

## Residual visualization

# Arrange data as layers (values & residuals)
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

    # Residual map subpanel
    res_map = plot_layer(
        layer;
        c=recentergrad(:PuOr, lims; rev=true),
        title="Residuals map",
        colorbar_title="Deviance residuals",
        clims=lims,
    )

    # Residual histogram subpanel
    res_hist = histogram(
        layer;
        bins=50,
        c=:PuOr,
        legend=:none,
        title="Residuals distribution",
        xlabel="Relative frequency",
        ylabel="Deviance residuals",
        orientation=:horizontal,
        ylims=lims,
        normalize=:probability
    )

    # Combine in single plot
    l = @layout [a{0.6w} b{0.38w}]
    res_plot = plot(
        res_map,
        res_hist;
        size=(850, 340),
        layout=l,
        bottommargin=3.0mm,
        dpi=200,
        rightmargin=[5.0mm 0mm],
        leftmargin=[5.0mm 5.0mm],
        kw...,
    )
    return res_plot
end

# Richness residuals plots
# richness_resplot residuals_plot(richres_layer; title = "Richness Poisson GLM")
richness_qp_resplot = residuals_plot(richres_qp_layer; title="Richness Quasipoisson GLM")
richness_nb_resplot = residuals_plot(
    richres_nb_layer;
    title="Richness Negative Binomial GLM",
    # yticks = [:auto :auto (-3:1:4, string.(-3:1:4))]
    yticks=[:auto (-3:1:4, string.(-3:1:4))],
)

# LCBD residuals plots
lcbd_resplot = residuals_plot(lcbdres_layer; title="LCBD Gamma GLM")
lcbd_br_resplot = residuals_plot(lcbdres_br_layer; title="LCBD Beta Regression")

# Combine richness & LCBD plots
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
xlims!(combined_resplot[2], (-0.003, 0.11))
xlims!(combined_resplot[4], (-0.003, 0.11))

# Save figures
if (@isdefined save_additional_figures) && save_additional_figures == true
    savefig(combined_resplot, joinpath("fig", "bart", "09_bart_residuals.png"))
end
