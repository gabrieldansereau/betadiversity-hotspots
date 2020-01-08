Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

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
mkdir coverfraction-csv
mkdir coverfraction/bare
mkdir coverfraction/crops
mkdir coverfraction/grass
mkdir coverfraction/moss
mkdir coverfraction/shrub
mkdir coverfraction/snow
mkdir coverfraction/tree
mkdir coverfraction/urban
mkdir coverfraction/water-permanent
mkdir coverfraction/water-seasonal

# Batch commands for each land cover variable
for i in $(ls coverfraction)
do
    # Copy landcover data in 1 folder per variable
    for j in landcover*; do cp "$j"/*"$i"-coverfraction-layer* coverfraction/"$i"/; done
    # Set resolution to 10 arc-minutes, keep 255 as no data value, merge all layers in one
    gdalwarp -tr 0.166666666666666650 0.1666666666666666570 -r average coverfraction/"$i"/W*.tif coverfraction/"$i"/"$i"-warp.tif
    # Transform layer to .asc format, which can be opened in text editor
    gdal_translate -of AAIGrid coverfraction/"$i"/"$i"-warp.tif coverfraction/"$i"/"$i".asc
    # Remove first 6 lines of .asc file, creates grid that can be loaded easily
    tail -n +7 coverfraction/"$i"/"$i".asc > ../BioClim/assets/landcover/lc_"$i"_10m.csv
done
=#

##  Test landcover variables

# Define coordinates range
lon_range = (-145.0, -50.0)
lat_range = (20.0, 75.0)

# Get worldclim variables
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);

# Test loading variables
lc_layers = load_landcover(lon_range, lat_range)

# Plot layers
plotSDM(lc_layers[1])
plotSDM(wc_vars[1])

# Test for sites with landcover over 100
nul_layer = zeros(Float64, size(lc_layers[1].grid))
for l in lc_layers
    global nul_layer += l.grid
end
nul_layer
nul_layer_nonan = filter(!isnan, nul_layer)
describe(nul_layer_nonan)
filter(x -> x > 110, nul_layer_nonan)
filter(x -> x > 120, nul_layer_nonan)
