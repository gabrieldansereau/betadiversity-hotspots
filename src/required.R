# Load packages
library(conflicted)
library(tidyverse)
library(here)
library(rgdal)
library(raster)
library(embarcadero)
library(viridis)
library(furrr)
library(ranger)
library(caret)
library(pbapply)
library(broom)
library(MASS)
library(betareg)
library(SpatialPack)

# Select parallel processing option
if (future::supportsMulticore()) {
    future::plan(future::multicore) # Preferred option when possible (Linux & Mac, not on RStudio)
} else {
    future::plan(future::multisession) # For Windows and RStudio (also works on Linux & Mac)
}

# Resolve conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("select", "dplyr")

# Load custom functions
source(here("src", "lib", "bart.R"))
