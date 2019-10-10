# BioClim SDM

This repository contains work for my M.Sc. on the identification of beta-diversity hotspots using species distribution models (SDMs). It is based on previous proof of concept by @tpoisot, my supervisor, at https://gitlab.com/tpoisot/BioClim.

This project is implemented in *Julia v1.2.0*. The required packages are listed in `src/required.jl`, and include the unreleased version of the `SimpleSDMLayers` package, available at https://github.com/EcoJulia/SimpleSDMLayers.jl.

The data used in this project comes from the *eBird Basic Dataset*. The project is for now focused on all **Warblers species** (*Parulidae* family) in North America (CA, US, MX)
```
eBird Basic Dataset. Version: EBD_relJun-2019. Cornell Lab of Ornithology, Ithaca, New York. Jun 2019.
```
Note however that neither the data nor the SDM prediction results are hosted in this remote repository, due to size limitations. Further updates will make everything accessible and reproducible once a proper solution is found.

The initial project used to download data through the `GBIF` package. Although some functions still support the `GBIFRecords` type, it has been abandonned in the analyses.

## Repository structure

The repository is organized has follows:

* `archive/` contains outdated elements kept in case of future need.

* `assets/` contains the `worldclim` data downloaded through the `SimpleSDMLayers` package.

* `data/` is used locally to store the data.
   * `jld2/` contains exported *Julia* `.jld2` elements, as the SDM predictions.
   * `proc/` contains processed CSV data.
   * `raw/` contains the raw CSV datasets from eBird.

* `fig/` contains the figures produced.
   * `raw/` and `sdm/` contain result figures, see analysis workflow for details.

* `pres/` contains presentations made about the project.

* `src/` contains all the scripts used in the project. Ordered scripts in this directory are the main represent the main steps of the analyses. Other general-use scripts are here as well. Subfolders contain scripts with a more specific use.
   * `lib/` is the library of the custom functions used.
   * `test/` contains random testing scripts kept in case of future use.
   * `raw/` and `sdm/` contain the ordered analysis scripts, see analysis workflow for details.

## Analysis workflow

The general workflow of the analyses is as follows:

1. `src/01_data-preparation.jl` is used to prepare the data in `data/raw` and saves the results in `data/proc`.

1. `src/02_master-raw.jl` is used to call the subscripts in `src/raw/` and perform the analyses on the raw observation data, without any SDM prediction. Results are stored in `data/jld2/` and `fig/raw/`.

1. `src/03_master-sdm.jl` is used to call the subscripts in `src/sdm/`, which makes the SDM predictions and performes the analyses on them. Results are stored in `data/jld2` and `fig/raw/`

The `master` scripts are meant to reproduce all the results, but do not overwrite figures everytime. I often prefer to work directly in the subscripts for that.

## Main results

### Single species

`src/raw/01_raw_presence-absence.jl` and `src/raw/01_sdm_predictions.jl` map the observed and predicted distribution of single species, respectively, across my range of interest.`src/sdm/sdm_single_species.jl` is used the predict the species' distibution across their whole range.

![Single species - raw][raw_single-sp] ![Single species - sdm][sdm_single-sp]

[raw_single-sp]: fig/raw/01_raw_sp-Setophaga_petechia.pdf
[sdm_single-sp]: fig/sdm/01_sdm_sp-Setophaga_petechia.pdf

### Species richness

`src/raw/03_raw_richness.jl` and `src/sdm/03_sdm_richness.jl` map the observed and predicted species richness (number of species per site).

![Species richness - raw][raw_richness] ![Species richness - sdm][sdm_richness]

[raw_richness]: fig/raw/03_raw_richness.pdf
[sdm_richness]: fig/sdm/03_sdm_richness.pdf

### LCBD

`src/raw/05_raw_lcbd.jl` and `src/sdm/05_sdm_lcbd.jl` map the observed and predicted LCBDs per site.

![LCBD - raw][raw_lcbd] ![LCBD - sdm][sdm_lcbd]

[raw_lcbd]: fig/raw/05_raw_lcbd.pdf
[sdm_lcbd]: fig/sdm/05_sdm_lcbd.pdf

### Richness-LCBD relationship

`src/raw/06_raw_relation-lcbd-richness.jl` and `src/sdm/06_sdm_relation-lcbd-richness.jl` plot the relationship between species richness and LCBD per site for observed and predicted data.

![relationship - raw][raw_relationship] ![relationship - sdm][sdm_relationship]

[raw_relationship]: fig/raw/06_raw_relation-lcbd-richness.pdf
[sdm_relationship]: fig/sdm/06_raw_relation-lcbd-richness.pdf

## Other useful scripts

* `02_xxx_Y-matrix.jl` is used to calculate and plot the site by species matrix Y.

* `04_xxx_community.jl` is used to calculate Pielou's diversity index.

## Details on library scripts

This code is built around the `SimpleSDMLayer` types, which are used to store the
worldclim variables, and also store the ouput of the prediction. Additional support
functions are in `src/lib/`.

1. `src/lib/SDMLayer.jl` has some utility functions to extend those from `SimpleSDMLayers`, notably to manipulate and integrate `DataFrames`.

1. `src/lib/shapefiles.jl` has a function to download the shapefiles for plotting,
and a function to clip them so that they overlap with a `SDMLayer`.

1. `src/lib/csvdata.jl` has functions to prepare data extracted from CSV files.

1. `src/lib/bioclim.jl` has the actual functions to perform the SDMs.

1. `src/lib/beta-div.jl` has functions to compute beta-diversity statitstics.

1. `src/lib/plotSDM.jl` contains a function to allow easier plotting of the `SimpleSDMLayer` type elements.
