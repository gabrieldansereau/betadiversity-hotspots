# BioClim SDM

This repository contains **proof of concept** code to implement the *bioclim*
species distribution model in Julia 1.0. It requires `GDAL`, `Plots`, and
the unreleased version of the `GBIF` package (everything you need is in the
`Project.toml` / `Manifest.toml` files).

## Overview

This code is built around the `SDMLayer` type, which is used to store the
worldclim variables, and also store the ouput of the prediction. The support
functions are in `lib/`.

1. `lib/SDMLayer.jl` has the type and some utility functions, notably to
manipulate and integrate `GBIFRecords`.

2. `lib/gdal.jl` has a function to read GeoTiff files. It is probably *not*
the most elegant way of reading them, but It Works Well Enough.

3. `lib/worldclim.jl` has functions to get the worldclim data from a folder,
and that's pretty much it.

4. `lib/bioclim.jl` is the actual function to perform the SDM. It has a
non-optimal function to find the quantiles, where most of the time is spent.

## Example code

See in `main.jl` -- this will get up to about 800 occurrences of the american
robib (*Turdus migratorius*) in Canada, and display the prediction. The output
is cropped so that only the values above the 10th percentile are returned.

![Example of the SDM][sdm]

[sdm]: sdm.png
