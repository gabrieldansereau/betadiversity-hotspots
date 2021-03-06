#!/bin/bash

## Bash commands to download & prepare landcover data from Zenodo
# Just run following in bash terminal:
# bash src/shell/landcover_coarsen.sh

# Make sure file is executable. If not (e.g. Permission denied error), run:
# chmod +x src/shell/landcover_coarsen.sh

####

## Go to landcover directory
cd assets/landcover/
landcover_variables=(bare crops grass moss shrub snow tree urban water-permanent water-seasonal)

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

## Fix layer extents
for i in "${landcover_variables[@]}"
do
    # Rename layers to allow overwrite
    mv lc_"$i"_10m.tif lc_"$i"_10m_previous.tif
    # Set layer extent to exact values
    gdalwarp -te -180.0 -60.0 180.0 80.0 lc_"$i"_10m_previous.tif lc_"$i"_10m.tif
    # Remove previous layers
    rm lc_"$i"_10m_previous.tif
done