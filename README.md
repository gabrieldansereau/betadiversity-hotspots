# BioClim SDM

This repository contains **proof of concept** code to implement the *bioclim*
species distribution model in Julia 1.0. It requires `GDAL`, `Shapefiles`,
`Plots`, and the unreleased version of the `GBIF` package (everything you
need is in the `Project.toml` / `Manifest.toml` files).

## Overview

This code is built around the `SDMLayer` type, which is used to store the
worldclim variables, and also store the ouput of the prediction. The support
functions are in `lib/`.

1. `lib/SDMLayer.jl` has the type and some utility functions, notably to
manipulate and integrate `GBIFRecords`.

1. `lib/gdal.jl` has a function to read GeoTiff files. It is probably *not*
the most elegant way of reading them, but It Works Well Enough.

1. `lib/shapefiles.jl` has a function to download the shapefiles for plotting,
and a function to clip them so that they overlap with a `SDMLayer`.

1. `lib/worldclim.jl` has functions to get the worldclim data from a folder,
and that's pretty much it.

1. `lib/bioclim.jl` is the actual function to perform the SDM.

## Example code

### Single species

See in `main.jl` -- this will get up to about 3000 occurrences of some species
in Canada and the US from [GBIF], and display the prediction. The output is
cropped so that only the values above the 10th percentile are returned. Most
of the code in `main.jl` is actually the downloading of data, and the plotting.

[GBIF]: http://gbif.org

The time required to run the SDM scales with the size of the grid (linearly
with the number of cells). On a standard laptop, using the 10m worldclim data,
generating a prediction for the entire US + Canada takes about one second.

![Example of the SDM][sdm]

[sdm]: sdm.png

### Multiple species

In `community.jl`, there is some code to generate SDMs for species of Warblers,
and then aggregate the output. The output is, specifically, Pielou's evenness
of every SDM score, so that higher values means equal habitat suitability,
and lower values means unequal habitat suitability across species. The
scores are corrected for the number of species *locally present*. *i.e.*
species with a non-`NaN` suitability at the given pixel.

![Example of the SDM on warblers][warblers]

[warblers]: warblers.png
