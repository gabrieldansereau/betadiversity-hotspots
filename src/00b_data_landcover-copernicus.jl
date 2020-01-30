import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Conditional arguments
# save_figures = true # optional

## Bash commands to download & prepare data, to run in terminal in ../landcover/
#=
cd ~/github/landcover/
# Download Copernicus land cover data, many downloads to cover whole extent
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W160N80_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W160N80.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W160N60_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W160N60.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W160N40_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W160N40.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W140N80_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W140N80.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W140N60_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W140N60.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W140N40_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W140N40.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W120N80_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W120N80.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W120N60_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W120N60.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W120N40_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W120N40.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W100N80_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W100N80.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W100N60_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W100N60.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W100N40_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W100N40.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W080N80_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W080N80.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W080N60_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W080N60.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W080N40_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W080N40.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W060N80_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W060N80.zip
wget https://s3-eu-west-1.amazonaws.com/vito-downloads/ZIPfiles/W060N60_ProbaV_LC100_epoch2015_global_v2.0.1_products_EPSG-4326.zip -O landcover_W060N60.zip

# Unzip files in separate directories
for i in *.zip; do unzip "$i" -d "${i%%.zip}"; done

# Delete zip files
rm *.zip

# Create repositories if needed
mkdir coverfraction
(cd coverfraction; mkdir bare crops grass moss shrub snow tree urban water-permanent water-seasonal)

# Batch commands for each land cover variable
for i in $(ls coverfraction/)
do
    # Copy landcover data in 1 folder per variable
    for j in landcover*; do cp "$j"/*"$i"-coverfraction-layer* coverfraction/"$i"/; done
    # Set resolution to 10 arc-minutes & merge all layers in one
    gdalwarp -tr 0.166667 0.166667 -r average coverfraction/"$i"/W*.tif ../BioClim/assets/landcover/lc_"$i"_10m.tif
    # Set resolution to 5 arc-minutes & merge all layers in one
    gdalwarp -tr 0.0833333 0.0833333 -r average coverfraction/"$i"/W*.tif ../BioClim/assets/landcover/lc_"$i"-5m.tif
done
=#

##  Test landcover variables

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
# save_figures = true
if (@isdefined save_figures) && save_figure == true
    savefig(wc_plot, "fig/00b_wc1-temperature.pdf")
    savefig(lc_plot, "fig/00b_lc2-crops.pdf")
    @info "Figures saved"
else
    @info "Figures not saved"
end
