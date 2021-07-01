# Avoid function names conflicts
library(conflicted)
# Load packages
library(betareg)
library(broom)
library(caret)
library(embarcadero)
library(furrr)
library(here)
library(MASS)
library(pbapply)
library(ranger)
library(raster)
library(rgdal)
library(SpatialPack)
library(tictoc)
library(tidyverse)
library(viridis)

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
