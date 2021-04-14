# Load packages
library(conflicted)
library(tidyverse)
library(here)
library(embarcadero)
library(viridis)
library(furrr)
library(ranger)
library(caret)
library(pbapply)
library(broom)
plan(multiprocess)

# Resolve conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("select", "dplyr")

# Load custom functions
source(here("src", "lib", "R", "bart.R"))