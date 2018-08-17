# BioClim SDM

This repository contains **proof of concept** code to implement the *bioclim*
species distribution model in Julia 1.0. It requires `GDAL`, `Plots`, and
the unreleased version of the `GBIF` package (everything you need is in the
`Project.toml` / `Manifest.toml` files).

This code is built around the `SDMLayer` type, which is used to store the
worldclim variables, and also store the ouput of the prediction. The support
functions are in `lib/`.

1. `lib/SDMLayer.jl` has the type and some utility functions, notably to
manipulate and integrate `GBIFRecords`.

2. 
