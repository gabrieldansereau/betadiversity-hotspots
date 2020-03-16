import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Conditional arguments
# save_figures = true # should figures be overwritten (optional)

## Run bash scripts to download & coarsen landcover data from Zenodo (if files missing)
lc_files = readdir("assets/landcover/")
# Check if landcover files are missing
if !any(startswith.(lc_files, r"^lc_"))
    # Check if full resolution files are missing
    if !any(startswith.(lc_files, r"^landcover_copernicus_global_100m"))
        # Download full resolution files
        # BEWARE, can be long to download, 25 GB of data
        run(`bash src/bin/landcover_download.sh`)
    end
    # Coarsen resolution
    run(`bash src/bin/landcover_coarsen.sh`)
end

## Test landcover variables
# Define coordinates range
lon_range = (-145.0, -50.0)
lat_range = (20.0, 75.0)

# Test loading variables
lc_vars = landcover(1:10, resolution = "10")
lc_vars = landcover(1:10, resolution = "5")
lc_vars = map(x -> landcover(x, resolution = "5")[lon_range, lat_range], 1:10)
fig1 = plotSDM(lc_vars[2])

# Plot worldclim to compare
@time wc_vars = pmap(x -> worldclim(x, resolution = "5")[lon_range, lat_range], 1:19);
fig2 = plotSDM(wc_vars[1])
fig1
fig2

# Test for sites with landcover over 100
nul_layer = zeros(Float64, size(lc_vars[1].grid))
for l in lc_vars
    global nul_layer += l.grid
end
nul_layer
nul_layer_nonan = filter(!isnan, nul_layer)
describe(nul_layer_nonan)
filter(x -> x > 110, nul_layer_nonan)
filter(x -> x > 120, nul_layer_nonan)

## Plot environmental variables examples
# Plot wcvars1 (temperature)
wc_plot = plotSDM(wc_vars[1], clim=extrema(filter(!isnan,wc_vars[1].grid)),
                  colorbar_title="Annual Mean Temperature (°C)", dpi=300)
# Plot lcvars2 (urban)
lc_plot = plotSDM(lc_vars[2], colorbar_title="Crops land cover (%)", dpi=300)

## Export figures
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved (environmental variables)"
    savefig(wc_plot, "fig/00b_wc1-temperature.png")
    savefig(lc_plot, "fig/00b_lc2-crops.png")
else
    @info "Figures not saved (environmental variables)"
end
