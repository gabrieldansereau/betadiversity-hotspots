# BioClim SDM

This repository contains work for my M.Sc. on the identification of beta-diversity hotspots using species distribution models (SDMs). It is based on previous proof of concept by @tpoisot, my advisor, at https://gitlab.com/tpoisot/BioClim.

This project is implemented in *Julia v1.3.1*. The required packages and versions are listed in `Project.toml`.

The data used in this project comes from the *eBird Basic Dataset*. The project is for now focused on all warblers species (*Parulidae* family) in North America (CA, US, MX).

```
eBird Basic Dataset. Version: EBD_relJun-2019. Cornell Lab of Ornithology, Ithaca, New York. Jun 2019.
```
Note however that neither the data nor the SDM prediction results are hosted in this remote repository, due to size limitations. Further updates will make everything accessible and reproducible once a proper solution is found.

The initial project used to download data through the `GBIF` package. Although some functions still support the `GBIFRecords` type, it has been abandoned in the analyses.

## Repository structure

The repository is organized has follows:

* `archive/` contains outdated elements kept in case of future need.

* `assets/` contains the _Worldclim 2.0_ climate data (downloaded through the `SimpleSDMLayers` package) and the _Copernicus_ landcover data (downloaded in `src/00b_data_landcover-copernicus.jl`)

* `data/` is used locally to store the data.
   * `jld2/` contains exported *Julia* `.jld2` elements, as the SDM predictions.
   * `proc/` contains processed CSV data.
   * `raw/` contains the raw CSV datasets from eBird.

* `docs/` contains documents such as reports and presentations about my project.

* `fig/` contains the figures produced.
   * `raw/` and `sdm/` contain result figures, see analysis workflow for details.

* `src/` contains all the scripts used in the project. Ordered scripts in this directory represent the main steps of the analyses. Subfolders contain scripts with a more specific use.
   * `lib/` is the library of the custom functions used.
   * `others/` contains useful scripts that are not part of the main analyses
   * `test/` contains random testing scripts kept in case of future use.

## Analysis workflow

All analysis scripts are in `src/`.

* `master.jl` can be used to run all the analyses and produce the figures.
* `required.jl` loads all the required packages and library functions.

Else, the general workflow of the analyses is as follows:

1. `src/00a_data_ebd-preparation.jl` is used to prepare the data from eBird in `data/raw` (pre-processed, script coming soon) and saves the results in `data/proc`.

1. `src/00b_data_landcover-copernicus.jl` is used to prepare the landcover data from Copernicus.

1. `src/01_distributions` is used to compute, map and save the species distributions, either for raw data or with SDM predictions, depending on the `outcome` variable.

1. `02_Y-matrix.jl` is used to transform the distributions into a community matrix Y (site x species) used in the following analyses.

1. `03_richness.jl` is used to calculate and map species richness per site.

1. `04_evenness.jl` is used to calculate Pielou's evenness index among sites.

1. `05_lcbd.jl` is used to calculate and map LCBD indices (local contributions to beta diversity) among sites.

1. `06_relationship_lcbd-richness.jl` is used to plot the relationship between sites species richness and LCBD indices.

## Main results

### Single species

![Single species - Raw][raw_single-sp]

![Single species - SDM][sdm_single-sp]

[raw_single-sp]: fig/raw/01_raw_sp-Setophaga_petechia.pdf
[sdm_single-sp]: fig/sdm/01_sdm_sp-Setophaga_petechia.pdf

### Species richness

![Species richness - Raw][raw_richness]

![Species richness - SDM][sdm_richness]

[raw_richness]: fig/raw/03_raw_richness.pdf
[sdm_richness]: fig/sdm/03_sdm_richness.pdf

### LCBD

![LCBD - Raw][raw_lcbd]
![LCBD - SDM][sdm_lcbd]

[raw_lcbd]: fig/raw/05_raw_lcbd-transf.pdf
[sdm_lcbd]: fig/sdm/05_sdm_lcbd.pdf

### Richness-LCBD relationship

![Relationship][relationship]

[relationship]: fig/06_relationship_lcbd-richness-transf.png

## Details on library scripts

This code is built around the `SimpleSDMLayer` types, which are used to store the
environmental variables, and also store the ouput of the prediction. Additional support
functions are in `src/lib/`, where:

1. `beta-div.jl` has functions to compute beta-diversity statitstics.

1. `bioclim.jl` has the actual functions to perform the BIOCLIM SDM model.

1. `csvdata.jl` has functions to prepare data extracted from CSV files.

1. `landcover.jl` has functions to extract and prepare the landcover data.

1. `overloads.jl` has some utility functions to extend those from `SimpleSDMLayers`, notably to manipulate and integrate `DataFrames`.

1. `plotSDM.jl` contains a function to allow easier plotting of the `SimpleSDMLayer` type elements.

1. `presence-absence.jl` has the function to convert the raw data into a presence-absence layer.

1. `shapefiles.jl` has a function to download the shapefiles for plotting,
and a function to clip them so that they overlap with a `SimpleSDMLayer`.
