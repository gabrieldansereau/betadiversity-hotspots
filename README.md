# Beta Diversity Hotspots

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6024392.svg)](https://doi.org/10.5281/zenodo.6024392)

### Evaluating ecological uniqueness over broad spatial extents using species distribution modelling

![Effect of scaling - BARTs][bart_scaling]

[bart_scaling]: fig/bart//05_bart_scaling.gif

This repository contains work for my *M.Sc.* on the identification of beta diversity hotspots using species distribution models (SDMs). The results are part of a manuscript available as a [preprint on EcoEvoRxiv](https://ecoevorxiv.org/tvmyg). There is also [a specific repository](https://github.com/gabrieldansereau/ms_betadiversity_hotspots) for the manuscript.

This project is implemented in _Julia v1.6.1_. The required packages and versions are listed in `Project.toml`. To install them, run the first lines of [src/required.jl](./src/required.jl). Some steps are also implemented in _R v4.1.0_, with packages & versions tracked by `renv`. More details below. 

The data used in this project comes from the [eBird Basic Dataset](https://ebird.org/science/use-ebird-data/download-ebird-data-products) from June 2019. The project is for now focused on all warblers species (*Parulidae* family) in North America (CA, US, MX).

    eBird Basic Dataset. Version: EBD_relJun-2019. Cornell Lab of Ornithology,  Ithaca, New York. Jun 2019.

Note however that the data is not hosted in this remote repository due to size limitations.

## Repository structure

The repository is organized as follows:

- `assets/` contains the pre-coarsened *Copernicus* land cover data (downloaded and coarsened in `src/00c_data_landcover-copernicus.jl`).

- `data/` is used to store the data.
  - `jld2/` contains exported *Julia* `.jld2` elements, such as SDM predictions. Earlier versions relied heavily on these, but since they were too large to be version controlled, they have now been replaced by raster files in `raster/`. Some `.jld2` are still exported here and there but are not central to the analyses.
  - `proc/` contains processed CSV data. Importantly, the prepared eBird data and some BART predictions are locally stored here but are not version controlled due to their size.
  - `raster/` contains raster files with species distributions (observed and predicted) and environmental data layers. Raster files are now central to the workflow and are used to save and reload data between scripts.
  - `raw/` contains the raw CSV datasets from eBird (not version-controlled).
  - `rdata/` contains `.RData` files used as backups in the R scripts, which are not essential and are not version controlled.

- `fig/` contains the figures produced, organized by outcome (`bart` for figures based on predicted data and `raw` for figures based on observed data).

- `src/` contains all the scripts used in the project. Ordered scripts in this directory represent the main steps of the analyses. Subfolders contain scripts with a more specific use.
  - `lib/` is the library of the custom functions used in the main scripts.
  - `others/` contains useful scripts that are not part of the main analyses
  - `shell/` contains Bash scripts used for some operations.

## Analysis workflow

All analysis scripts are in `src/`.

- `main.jl` can be used to run all the analyses and produce the figures.
- `required.jl` loads all the required packages and library functions.

Else, the general workflow of the analyses is as follows:

1. `00a_ebd_extraction.jl` extracts the Warblers data from the complete EBD to `data/raw` (not version controlled).

1. `00b_ebd_preparation.jl` prepares the Warblers data in `data/raw` for the analyses, then saves the results in `data/proc` (not version-controlled)

1. `00c_landcover.jl` prepares the landcover data from Copernicus and exports the environmental data as CSVs in `data/proc` and as TIFF files in `data/raster`.

1. `01_distributions` assembles the species distributions from the raw data as layers, then exports these to `data/raster`. It also produces examples of single species maps.

1. `02_training_bart.R` trains BARTs (Bayesian Additive Regression Trees) in _R_ (package `embarcadero`) based on the distribution and environmental rasters, then predicts the species distributions (exported as CSV files, which are not version controlled).

1. `03_predictions_bart.jl` assembles the predicted distributions as layers and exports them as raster files.

1. `04_full-extent.jl` performs the main analysis steps (on the full spatial extent): getting species richness and LCBD values per site and verifying the relationship between the two. These steps can be performed on either the observed or predicted distributions.

1. `05_subareas.jl`  reapplies the analyses on smaller regions and investigates the effect of the spatial scale on the results.

1. `06_moving-windows.jl` investigates the effect of the proportion of rare species on the relationship between species richness and LCBD values at varying scales.

1. `07_comparison_data.jl` re-runs the main analysis steps on both the observed and predicted data and prepares the results for comparison in the following scripts.

1. `08_comparison_glm.jl` performs GLMs in _R_ to compare the observed and predicted results and saves the results to be plotted in the next script.

1. `09_comparison_plots.jl` produces plots comparing the observed and predicted results. The comparison is made by comparing the results directly (called difference plots) or the GLM residuals produced in the previous script (called residual plots).

## Main results

### Observed vs predicted values

![Comparison of observed & predicted results][comparison_results]

[comparison_results]: fig/bart/09_bart_combined.png

### Subareas

![Subareas - BARTs][bart_subareas]

[bart_subareas]: fig/bart/05_bart_subareas.png

### Effect of scaling

![Effect of scaling - BARTs][bart_scaling]

[bart_scaling]: fig/bart//05_bart_scaling.gif

## Details on library scripts

This code is built around the package [SimpleSDMLayers.jl](https://github.com/EcoJulia/SimpleSDMLayers.jl) and its `SimpleSDMLayer` types, which are used to store the environmental variables and the species distributions.

1. `analysis.jl` contains the functions to perform the main analyses.

1. `bart.R` contains utility functions for the BART analyses.

1. `betadiv.jl` contains functions to compute beta diversity statistics.

1. `csvdata.jl` contains functions to prepare the data extracted from CSV files.

1. `landcover.jl` contains functions to extract and prepare the landcover data (similarly to the other data sources in `SimpleSDMLayers.jl`).

1. `plotting.jl` contains a function to allow easier plotting of the `SimpleSDMLayer` type elements.

1. `presence-absence.jl` contains the function to convert the raw data into a presence-absence layer.

1. `shapefiles.jl` contains a function to download the background shapefiles for plotting, and a function to clip them so that they overlap with a `SimpleSDMLayer`.

2. `version-control.jl` contains the list of important files which are too large to be version controlled and a set of custom functions to track their changes. See additional notes for details.

## Additional notes

- For each important file that is too large to be version controlled, a version-controlled placeholder file was created (for example `data/proc/ebd_warblers_prep_placeholder.csv`) to record the time where the large file was last updated on purpose. The placeholder is updated when the files are changed, and the functions will trigger a warning prompting to make sure the change was made on purpose. If it was, the placeholder should be re-committed with the new modification time. If the file was overwritten without a change (and the user is sure of it), the placeholder change can be discarded.
- This project is based on previous proof of concept by [@tpoisot](https://github.com/tpoisot), my M.Sc. advisor, at <https://gitlab.com/tpoisot/BioClim>.
