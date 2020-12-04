# Load data
outcome = "bart"
outcome = "raw"
include(joinpath(pwd(), "src", "04_analysis.jl"))

# Load temperature
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
temp = worldclim(1)[coords]
# temp5 = worldclim(1, resolution = 5.0)[coords]

# Clip to QC
coords_qc = (left = -82.5, right = -51.0,
             bottom = 41.0, top = 71.0)
temp_qc = temp[coords_qc]
plot(temp_qc)
# Set Greenland values to nothing
size(temp_qc)
xmin = SimpleSDMLayers._match_longitude(temp_qc, -55.0)
xmax = SimpleSDMLayers._match_longitude(temp_qc, coords_qc.right)
ymin = SimpleSDMLayers._match_latitude(temp_qc, 60.0)
ymax = SimpleSDMLayers._match_latitude(temp_qc, coords_qc.top)
temp_qc.grid[ymin:ymax, xmin:xmax] .= nothing
plot(temp_qc)
# Default to QC values
temp_full = copy(temp)
temp = temp_qc

# Temperature map
temp_plot = plot(temp_qc)

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
removebackground!(temp_plot)
# Save results
savefig(plot(temp_plot, dpi = 200), "./docs/qcbs/2020/fig/qc_none_inferno.png")

## Paint it white
# White inside background
temp_white = removebackground!(plot(temp_qc), bg_inside = :white)
# Save results
savefig(plot(temp_white, dpi = 200), "./docs/qcbs/2020/fig/qc_inside_inferno.png")

## Paint it viridis
# Viridis colorpalette
temp_viridis = removebackground!(plot(temp_qc, c = :viridis))
savefig(plot(temp_viridis, dpi = 200), "./docs/qcbs/2020/fig/qc_none_viridis.png")
# White inside
temp_viridis_white = removebackground!(plot(temp_qc, c = :viridis), bg_inside = :white)
savefig(plot(temp_viridis_white, dpi = 200), "./docs/qcbs/2020/fig/qc_inside_viridis.png")

## Create empty background map
empty_plot = plot(temp_qc, c = :lightgrey)
removebackground!(empty_plot)
savefig(plot(empty_plot, dpi = 200), "./docs/qcbs/2020/fig/qc_none_lightgrey.png")
savefig(plot(empty_plot, dpi = 200, bg_inside = :white), "./docs/qcbs/2020/fig/qc_inside_lightgrey.png")
# No background white
empty_plot_white = removebackground!(plot(temp_qc, c = :white))
savefig(plot(empty_plot_white, dpi = 200), "./docs/qcbs/2020/fig/qc_none_white.png")


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
savefig(plot(rectplot, dpi = 200), "./docs/qcbs/2020/fig/subarea-map.png")
