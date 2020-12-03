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
