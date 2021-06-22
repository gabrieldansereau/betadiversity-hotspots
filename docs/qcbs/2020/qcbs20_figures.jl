## Prepare layers & data
import Pkg
Pkg.activate(".")
@time include(joinpath("..", "..", "..", "src", "required.jl")) # best way for Intellisense
figdir = joinpath("docs", "qcbs", "2020", "fig")

# Load temperature
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
temp_full = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
# temp5 = SimpleSDMPredictor(WorldClim, BioClim, 1; resolution = 5.0, coords...)

# Clip to QC
coords_qc = (left = -82.5, right = -51.0,
             bottom = 41.0, top = 71.0)
temp_qc = temp_full[coords_qc]
# Verify dimensions
plot(temp_qc)
# Set Greenland values to nothing
xmin = SimpleSDMLayers._match_longitude(temp_qc, -55.0)
xmax = SimpleSDMLayers._match_longitude(temp_qc, coords_qc.right)
ymin = SimpleSDMLayers._match_latitude(temp_qc, 60.0)
ymax = SimpleSDMLayers._match_latitude(temp_qc, coords_qc.top)
temp_qc.grid[ymin:ymax, xmin:xmax] .= nothing
plot(temp_qc)

## Create maps for QC
qc_inferno = plot(temp_qc)
qc_lightgrey = plot(temp_qc, c = :lightgrey)
qc_viridis = plot(temp_qc, c = :viridis)
qc_white = plot(temp_qc, c = :white)

## Remove background
# Custom function
function removebackground!(p; kw...)
    plot!(p,
          frame = :origin,
          colorbar = :none,
          guides = "",
          grid = :none,
          bg = :transparent,
          bg_inside = :transparent;
          kw...
          )
    return p
end

# Remove background
removebackground!(qc_inferno)
removebackground!(qc_lightgrey)
removebackground!(qc_viridis)
removebackground!(qc_white)

# Save results
savefig(plot(qc_inferno,   dpi = 200), joinpath(figdir, "qc_no-bg_inferno.png"))
savefig(plot(qc_lightgrey, dpi = 200), joinpath(figdir, "qc_no-bg_lightgrey.png"))
savefig(plot(qc_viridis,   dpi = 200), joinpath(figdir, "qc_no-bg_viridis.png"))
savefig(plot(qc_white,     dpi = 200), joinpath(figdir, "qc_no-bg_white.png"))

# Paint it white
plot!(qc_inferno, bg_inside = :white)
plot!(qc_lightgrey, bg_inside = :white)
plot!(qc_viridis, bg_inside = :white)

# Save results
savefig(plot(qc_inferno,   dpi = 200), joinpath(figdir, "qc_inside-bg_inferno.png"))
savefig(plot(qc_lightgrey, dpi = 200), joinpath(figdir, "qc_inside-bg_lightgrey.png"))
savefig(plot(qc_viridis,   dpi = 200), joinpath(figdir, "qc_inside-bg_viridis.png"))

## Draw rectangle over subareas
# Get coordinates
coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
coords_SW = (left = -120.0, right = -100.0, bottom = 30.0, top = 40.0)

# Rectangle coordinates
function rectangle_from_coords(xb,yb,xt,yt)
    [
        xb  yb
        xt  yb
        xt  yt
        xb  yt
        xb  yb
        NaN NaN
    ]
end
rect_NE = rectangle_from_coords(coords_NE.left,  coords_NE.bottom,
                                coords_NE.right, coords_NE.top)
rect_SW = rectangle_from_coords(coords_SW.left,  coords_SW.bottom,
                                coords_SW.right, coords_SW.top)

# Plot subarea rectangles
rectplot = plot(temp_full, c = :lightgrey, cb = :none)
plot!(rectplot, rect_NE[:, 1], rect_NE[:, 2], label = "NE subarea")
plot!(rectplot, rect_SW[:, 1], rect_SW[:, 2], label = "SW subarea")
plot!(rectplot, bg_outside = :transparent)
savefig(plot(rectplot, dpi = 200), joinpath(figdir, "subarea_map.png"))

## Remove background on analysis figures
# Predefine desired outcome
outcome = "bart"
# Load analysis script
include(joinpath("..", "..", "..", "src", "04_full-extent.jl"))

# Remove backgrounds
plot!(richness_plot, bg_outside = :transparent, title = "")
plot!(lcbdtr_plot, bg_outside = :transparent, title = "")
plot!(rel2d_plot, bg_outside = :transparent, title = "")

# Save results
savefig(plot(richness_plot, dpi = 200), joinpath(figdir, "res_richness.png"))
savefig(plot(lcbdtr_plot, dpi = 200), joinpath(figdir, "res_lcbd.png"))
savefig(plot(rel2d_plot, dpi = 200), joinpath(figdir, "res_relationship.png"))

## Distribution figures
# Single species distributions
spplot1 = plotSDM2(distributions[1][coords_qc], c = :BuPu)
spplot2 = plotSDM2(distributions[17][coords_qc], c = :BuPu)
spplot3 = plotSDM2(distributions[22][coords_qc], c = :BuPu)
# Remove backgrounds
removebackground!(spplot1)
removebackground!(spplot2)
removebackground!(spplot3)
# Save
savefig(plot(spplot1, dpi = 200), joinpath(figdir, "sp_bart_1.png"))
savefig(plot(spplot2, dpi = 200), joinpath(figdir, "sp_bart_2.png"))
savefig(plot(spplot3, dpi = 200), joinpath(figdir, "sp_bart_3.png"))

## Raw distribution figures
# Load data
distributions_bart = copy(distributions)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
distributions_raw = copy(distributions)
distributions = distributions_bart
# Plot
spplot1_raw = plotSDM2(distributions_raw[1][coords_qc], c = :BuPu)
spplot2_raw = plotSDM2(distributions_raw[17][coords_qc], c = :BuPu)
spplot3_raw = plotSDM2(distributions_raw[22][coords_qc], c = :BuPu)
# Remove backgrounds
removebackground!(spplot1_raw)
removebackground!(spplot2_raw)
removebackground!(spplot3_raw)
# Save
savefig(plot(spplot1_raw, dpi = 200), joinpath(figdir, "sp_raw_1.png"))
savefig(plot(spplot2_raw, dpi = 200), joinpath(figdir, "sp_raw_2.png"))
savefig(plot(spplot3_raw, dpi = 200), joinpath(figdir, "sp_raw_3.png"))

## Subarea comparison figure
include(joinpath("..", "..", "..", "src", "05_subareas.jl"))
plot!(combined_plot, bg_outside = :transparent)
savefig(plot(combined_plot, dpi = 200), joinpath(figdir, "subarea_comparison.png"))

## Relationship plot
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
outcome = "raw"
include(joinpath("..", "..", "..", "src", "04_full-extent.jl"))
plot(richness, lcbd, colour = :transparent, smooth = true)

# Arrange data
resdf = DataFrame([richness, lcbd])
rename!(resdf, :x1 => :richness, :x2 => :lcbd)
filter!(x -> !isnothing(x.richness) && !isnothing(x.lcbd), resdf)
resdf.richness = Array{Float64}(resdf.richness)
resdf.lcbd = Array{Float64}(resdf.lcbd)

# Fit loess model
using Loess
function loess_vals(x, y; kw...)
    model = loess(x, y; kw...)
    xloess = range(extrema(x)...; step = 0.1)
    yloess = Loess.predict(model, xloess)
    return xloess, yloess
end
# Plot
xs, ys = loess_vals(resdf.richness, resdf.lcbd)
plot(xs, ys)

# Blank relationship plot
x = 0:0.1:10.0
y = @__dot__ 0.001*x^2 - 0.041*x - 3.394 # Heino 2017
y = @__dot__ 0.002*x^2 - 0.041*x - 3.394 # Heino 2017
y = @__dot__ 0.025*x^2 - 0.442*x - 2.690 # Silva 2018
relplot = plot(x, y,
               c = :black,
               legend = false,
               ylabel = "LCBD",
               xlabel = "Richness",
               )
plot!(labelfontsize = 16,
      ylim = (-3.62, Inf),
      ticks = false,
      thickness_scaling = 2)
# plot!(xlim = (-0.25, 10.0),
#       ylim = (-4.75, Inf),
#       aspect_ratio = 4.5
#       ) # for Silva 2018
savefig(plot(relplot, dpi = 50, bg = :transparent),
        joinpath(figdir, "relationship-plot.png"))
