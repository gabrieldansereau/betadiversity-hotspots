import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Conditional arguments
# save_figures = true # should figures be overwritten (optional)

## Bash commands to download & prepare data, to run in terminal in ../landcover/
#=
cd ../landcover/
## Download Copernicus global land cover data from Zenodo
# BEWARE, can be very long, 25 GB of data in total
# Launching downloads in parallel, 1 core/variable = 10 cores at most, not much RAM needed
landcover_variables=(bare crops grass moss shrub snow tree urban water-permanent water-seasonal)
for i in "${landcover_variables[@]}"
do
    wget https://zenodo.org/record/3243509/files/ProbaV_LC100_epoch2015_global_v2.0.2_"$i"-coverfraction-layer_EPSG-4326.tif -O landcover_copernicus_global_100m_v2.0.2_"$i".tif &
done
wait
echo "Downloads - All done"
rm wget-log*

## Coarsen resolution
# BEWARE, can be very long and will take 10 cores, took 30 min and ~ 16GB of RAM in my case
for i in "${landcover_variables[@]}"
do
    # Set resolution to 10 arc-minutes
    gdalwarp -tr 0.166667 0.166667 -r average --config GDAL_CACHEMAX 500 -wm 500 -multi landcover_copernicus_global_100m_v2.0.2_"$i".tif lc_"$i"_10m.tif &
done
wait
echo "10 arc-minutes - All done"

for i in "${landcover_variables[@]}"
do
    # Set resolution to 5 arc-minutes
    gdalwarp -tr 0.0833333 0.0833333 -r average --config GDAL_CACHEMAX 500 -wm 500 -multi landcover_copernicus_global_100m_v2.0.2_"$i".tif lc_"$i"_5m.tif &
done
wait
echo "5 arc-minutes - All done"
=#

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
wc_plot = plotSDM(wc_vars[1], c=:auto)
heatmap!(wc_plot, clim=extrema(filter(!isnan,wc_vars[1].grid)),
         colorbar_title="Annual Mean Temperature (°C)", dpi=300)
# Plot lcvars2 (urban)
lc_plot = plotSDM(lc_vars[2], c=:auto)
heatmap!(lc_plot, colorbar_title="Crops land cover (%)", dpi=300)

## Export figures
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved (environmental variables)"
    savefig(wc_plot, "fig/00b_wc1-temperature.pdf")
    savefig(lc_plot, "fig/00b_lc2-crops.pdf")
else
    @info "Figures not saved (environmental variables)"
end
