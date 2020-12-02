outcome = "bart"
outcome = "raw"
include(joinpath(pwd(), "src", "04_analysis.jl"))

coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
temp = worldclim(1)[coords]
temp5 = worldclim(1, resolution = 5.0)[coords]

testplot = plotSDM2(lcbd, c = :viridis)
testplot = plotSDM2(richness, c = :viridis)
testplot = plotSDM2(temp)
testplot = plotSDM2(temp5)
testplot = plot(temp5)
testplot = plot(temp, c = :lightgrey) |> x -> plot!(richness, c = :viridis, clim = extrema(richness))
testplot = plot(lcbd, c = :viridis)

plot!(testplot,
      bg_outside = :transparent,
      # fg = :white
      )

plot!(
      frame = :origin,
      grid = :none,
      colorbar = :none,
      guides = "",
      bg = :transparent,
      )

savefig(plot(testplot, dpi = 200), "./fig/testplot.png")
savefig(plot(testplot, dpi = 200), "./fig/testplot2.png")
