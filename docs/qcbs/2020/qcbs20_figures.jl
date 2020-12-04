## Prepare layers & data
import Pkg
Pkg.activate(".")
@time include(joinpath("..", "..", "..", "src", "required.jl")) # best way for Intellisense

# Load temperature
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
temp_full = worldclim(1)[coords]
# temp5 = worldclim(1, resolution = 5.0)[coords]

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
savefig(plot(qc_inferno,   dpi = 200), "./docs/qcbs/2020/fig/qc_no-bg_inferno.png")
savefig(plot(qc_lightgrey, dpi = 200), "./docs/qcbs/2020/fig/qc_no-bg_lightgrey.png")
savefig(plot(qc_viridis,   dpi = 200), "./docs/qcbs/2020/fig/qc_no-bg_viridis.png")
savefig(plot(qc_white,     dpi = 200), "./docs/qcbs/2020/fig/qc_no-bg_white.png")

# Paint it white
plot!(qc_inferno, bg_inside = :white)
plot!(qc_lightgrey, bg_inside = :white)
plot!(qc_viridis, bg_inside = :white)

# Save results
savefig(plot(qc_inferno,   dpi = 200), "./docs/qcbs/2020/fig/qc_inside-bg_inferno.png")
savefig(plot(qc_lightgrey, dpi = 200), "./docs/qcbs/2020/fig/qc_inside-bg_lightgrey.png")
savefig(plot(qc_viridis,   dpi = 200), "./docs/qcbs/2020/fig/qc_inside-bg_viridis.png")

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
