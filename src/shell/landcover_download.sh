#!/bin/bash

## Bash commands to download & prepare landcover data from Zenodo
# Just run following in bash terminal:
# bash src/shell/landcover_download.sh

# Make sure file is executable. If not (e.g. Permission denied error), run:
# chmod +x src/shell/landcover_download.sh

####

## Go to landcover directory
cd assets/landcover/

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
