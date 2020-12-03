# Load data
outcome = "bart"
outcome = "raw"
include(joinpath(pwd(), "src", "04_analysis.jl"))

# Load temperature
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
temp = worldclim(1)[coords]
temp5 = worldclim(1, resolution = 5.0)[coords]

# Testing plots
testplot = plotSDM2(lcbd, c = :viridis)
testplot = plotSDM2(richness, c = :viridis)
testplot = plotSDM2(temp)
testplot = plotSDM2(temp5)
testplot = plot(temp5)
testplot = plot(temp, c = :lightgrey) |> x ->
    plot!(richness, c = :viridis, clim = extrema(richness))
testplot = plot(lcbd, c = :viridis)
testplot = plotSDM2(temp, c = :lightgrey)

## Remove background
# Custom function
function removebackground!(p; kw...)
    plot!(p,
          frame = :origin,
          colorbar = :none,
          guides = "",
          grid = :none,
          bg = :transparent;
          kw...
          )
    return p
end

# Remove background
removebackground!(testplot)

# Save result
savefig(plot(testplot, dpi = 200), "./docs/qcbs/2020/fig/testplot.png")

## Keep inside background white
testplot2 = plot(testplot,
                 bg_inside = :white,
                 # bg_outside = :transparent,
                 # fg = :white
                 )
# Save result
savefig(plot(testplot2, dpi = 200), "./docs/qcbs/2020/fig/testplot2.png")

## Create background map
bgmap = plot(temp5, c = :lightgrey)
removebackground!(bgmap)
savefig(plot(bgmap, dpi = 200), "./docs/qcbs/2020/fig/background-map.png")

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
rectplot = plot(temp, c = :lightgrey, cb = :none)
plot!(rectplot, rect_NE[:, 1], rect_NE[:, 2], label = "NE subarea")
plot!(rectplot, rect_SW[:, 1], rect_SW[:, 2], label = "SW subarea")
savefig(plot(rectplot, dpi = 200), "./docs/qcbs/2020/fig/subarea-map.png")
