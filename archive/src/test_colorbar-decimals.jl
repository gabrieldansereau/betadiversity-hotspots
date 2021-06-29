include("../required.jl")

outcome = "bart"

## Prepare data

# Load presence-absence data for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
# Create matrix Y
Y = Ymatrix(distributions; transform=false)
# Create matrix of observed sites only
Yobs = _Yobs(Y)
# Get richness
richness = calculate_richness(Y, distributions[1])
# Get LCBD values
lcbd_rel = calculate_lcbd(Y, distributions[1]; transform=true, relative=true)
# Get non-relative values
lcbd_abs = calculate_lcbd(Y, distributions[1]; transform=true, relative=false)
round.(Float64.(extrema(lcbd_abs)); sigdigits=4)
lcbd = lcbd_abs

# Create a simpler MWE
mwe = rand(10, 10) / 2e5

## Option 1: Default plots
gr()
heatmap(mwe; c=:viridis, colorbar_title="LCBD value")
heatmap(lcbd; c=:viridis, colorbar_title="LCBD value")
plot_layer(lcbd; c=:viridis, colorbar_title="LCBD value")

#= Problems:
    - Colorbar title overlaps with colorbar ticks
   Pros:
    - That's what I have already
=#

## Option 2: Add scale in colorbar title
# Suggestion from https://discourse.julialang.org/t/plots-jl-clean-scientific-formatting/38112/2
gr()
heatmap(mwe * 1e5; c=:viridis, colorbar_title="LCBD value (x 10⁻⁵)")
layer1e5 = replace(similar(lcbd), 0.0 => 1e5)
heatmap(lcbd * layer1e5; c=:viridis, colorbar_title="LCBD value (x 10⁻⁵)")
plot_layer(lcbd * layer1e5; c=:viridis, colorbar_title="LCBD value (x 10⁻⁵)")

# Alternative to rescale
lcbd_scaled = lcbd_abs .* 1e5
extrema(lcbd_scaled)
lcbd_expo = rescale(lcbd_abs, (extrema(lcbd_scaled)))
collect(lcbd_expo)
isapprox(collect(lcbd_expo), lcbd_scaled)

#= Problems:
    - Still a bit of overlap, but could be ok
    - Not sure if it's considered a valid presentation form
   Pros:
    - Easier to fit with what I already have
=#

## Option 2: PyPlot
pyplot()
heatmap(mwe; c=:viridis, colorbar_title="LCBD value")
heatmap(lcbd; c=:viridis, colorbar_title="LCBD value")
heatmap(replace(lcbd.grid, nothing => NaN); c=:viridis, colorbar_title="LCBD value")
plot_layer(lcbd; c=:viridis, colorbar_title="LCBD value")

# Customize colorbar ticks
# New features from https://github.com/JuliaPlots/Plots.jl/pull/3346
# But only for PyPlot for now...
# GR extension discussed in https://github.com/JuliaPlots/Plots.jl/issues/3174
heatmap(rand(10, 10); colorbar_ticks=([0.3, 0.6], ["0.300", "0.60"]))
# How to check these possible attributes
plotattr(:Subplot)
plotattr("colorbar")

#= Problems:
    - super slow for layers (long to display)
    - looks terrible for layers
    - some weird lines on the heatmap, like a bad printer...
    - doesn't work with plot_layer
    - python
   Pros:
    - has scientific notation out-of-the-box
    - can customize colorbar ticks
=#

## Option 3: PlotlyJS (or Plotly)
plotlyjs()
heatmap(mwe; c=:viridis, colorbar_title="LCBD value")
heatmap(lcbd; c=:viridis, colorbar_title="LCBD value")
heatmap(lcbd; c=:viridis, colorbar_title="LCBD value", leftmargin=100px)
plot_layer(lcbd; c=:viridis, colorbar_title="LCBD value")
savefig("test-plotlyjs.png")
savefig("test-plotlyjs.pdf")

#= Problems:
    - micro sign, not sure how to remove
    - space between plot & colorbar (because of aspect_ratio)
        - can improve with leftmargin
    - terrible png quality, but pdf is OK and not too big
   Pros:
    - overall fast & still nice
    - works with plot_layer
=#

## Option 4: PlotlyJS.jl (not the same as the plotlyjs backend)
# From https://discourse.julialang.org/t/scientific-notation-in-plots-jl-colorbar/57793/9
using PlotlyJS: PlotlyJS
trace = PlotlyJS.heatmap(;
    z=replace(lcbd.grid, nothing => NaN)',
    colorscale="Viridis",
    colorbar_title="LCBD value",
    colorbar_exponentformat="power",
)
PlotlyJS.plot(trace)

#= Problems:
    - doesn't work with plot_layer
    - works differently from Plots, so I don't know how to tweak it
    - not sure if it can be used to combine plots either, which I need for subplots
   Pros:
    - scientific notation
=#

## Option 5: CairoMakie
using CairoMakie: CairoMakie
CairoMakie.activate!()
CairoMakie.heatmap(
    replace(lcbd.grid, nothing => NaN)'; c=:viridis, colorbar_title="LCBD value"
)

#= Problems:
    - don't really know how it works
    - no scientific notation
   Pros:
    - scientific notation
=#

## Option 6: Play with layouts
# https://discourse.julialang.org/t/single-colorbar-for-heatmap-grid-layout/50609/3
# https://discourse.julialang.org/t/plots-jl-shared-colorbar-with-subplots/47269/4
gr()
begin
    p1 = plot(lcbd; c=:viridis)
    p2 = plot(; frame=:none)
    annotate!(p2, 0.5, 0.5, text("LCBD value", 11, :center, 90.0))
    l = @layout [a b{0.01w}]
    plot(p1, p2; layout=l, leftmargin=[0.0mm -7.0mm])
end
begin
    p1 = plot_layer(lcbd; c=:viridis)
    p2 = plot(; frame=:none)
    annotate!(p2, 0.5, 0.5, text("LCBD value", 11, :center, 90.0))
    l = @layout [a b{0.01w}]
    plot(p1, p2; layout=l, leftmargin=[0.0mm -7.0mm])
end
# Comparison with relative values
begin
    p1 = plot_layer(lcbd_rel; c=:viridis)
    p2 = plot(; frame=:none)
    annotate!(p2, 0.5, 0.5, text("LCBD value", 11, :center, 90.0))
    l = @layout [a b{0.01w}]
    plot(p1, p2; layout=l, leftmargin=[0.0mm -20.0mm])
end
plot_layer(lcbd_rel; c=:viridis, colorbar_title="LCBD value")
# Just need to play with margins

#= Problems:
    - cannot format colorbar to similar units
   Pros:
    - Should work with everything I have if I simply define a new plot_layer function
=#
