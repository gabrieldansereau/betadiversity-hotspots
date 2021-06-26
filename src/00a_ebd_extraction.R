#### Explore eBird Dataset ####
# Mostly following tutorial from http://strimas.com/ebird-best-practices/

library(conflicted)
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(here)

# resolve namespace conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# Set ebd path
# auk::auk_set_ebd_path("~/github/data/ebd/raw/")
# Restart session

#### 1. Explore EBD Sample ####

# Load data
ebd_sample <- read_tsv(here("..", "data", "ebd", "raw", "ebd_sample.txt"))

# Fix names format
names(ebd_sample) = gsub(" ", ".", names(ebd_sample))
names(ebd_sample)

# Check data
head(ebd_sample)
table(ebd_sample$PROTOCOL.TYPE)
table(ebd_sample$SCIENTIFIC.NAME)

# Check taxonomy
names(ebird_taxonomy)
head(ebird_taxonomy)

# Extract Parulidae taxonomy information
parulidae_taxonomy <- filter(ebird_taxonomy, family == "Parulidae")
head(parulidae_taxonomy)
View(parulidae_taxonomy)

# Extract Parulidae species names (species only, no hybrid-issf-spuh-...)
parulidae_species <- parulidae_taxonomy %>%
  filter(category == "species") %>%
  pull(scientific_name)
parulidae_species

#### 2. Test Data Extraction with EBD Sample ####

# Get ebd file
ebd <- auk_ebd("ebd_sample.txt")

# Apply filters
ebd_filters <- ebd %>%
  auk_species(parulidae_species[1:97]) %>%
  auk_country(c("CA", "US", "MX")) %>%
  auk_complete()
ebd_filters
# Not working for more than 97 species at the time...
# Generates "Error running AWK command"
# [1:97], [98:110], [14:110] all work

# Export filtered data
f_ebd <- here("..", "data", "ebd", "processed", "ebd_sample_test.csv")
auk_filter(ebd_filters, file = f_ebd, overwrite = T)

# Test exported file
test <- read_tsv(f_ebd)
names(test) = gsub(" ", ".", names(test))
head(test)
table(test$SCIENTIFIC.NAME)
table(test$PROTOCOL.TYPE)

#### 3. Extract data from complete database ####

# Check worldwide distribution for a few species (manually from eBird website)
amN <- c("Myioborus miniatus", "Myioborus pictus",
     "Cardellina rubrifrons")
amS <- c("Myioborus albifrons", "Myioborus melanocephalus",
     "Myioborus ornatus", "Myioborus flavivertex",
     "Myioborus pariae", "Myioborus castaneocapilla",
     "Myioborus brunniceps", "Myiothlypis rivularis",
     "Myiothlypis nigrocristata", "Myiothlypis leucoblephara",
     "Myiothlypis luteoviridis")
amC <- c("Myioborus torquatus", "Cardellina versicolor",
     "Cardellina rubra")
nul <- c("Myioborus albifacies", "Myioborus cardonai")
# Remove a few species not present in North America (to have 97 at most and run command onloy once)
to_remove <- c(amS, nul)
fewer_parulidae_species <- parulidae_species[!(parulidae_species %in% to_remove)]
fewer_parulidae_species

# Get ebd file
ebd <- auk_ebd("ebd_relJun-2019.txt",
         file_sampling = "ebd_sampling_relJun-2019.txt")
ebd_sampling <- auk_sampling("ebd_sampling_relJun-2019.txt")
# Apply filters
ebd_filters <- ebd %>%
  auk_species(fewer_parulidae_species) %>%
  auk_country(c("CA", "US", "MX")) %>%
  auk_complete()
ebd_filters
ebd_sampling_filters <- ebd_sampling %>%
  auk_country(c("CA", "US", "MX")) %>%
  auk_complete()
ebd_sampling_filters
# Export filtered data !!! SEVERAL HOURS !!!!
f_ebd <- here("..", "data", "ebd", "processed", "ebd_warblers.csv")
f_sampling <- here("..", "data", "ebd", "processed", "ebd_warblers_sampling.csv")
auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
auk_filter(ebd_sampling_filters, file = f_sampling)
# File is too big, needs to cut in terminal

#### 4. Cut file using bash command ####

# Select variables to keep
vars <- c("GLOBAL.UNIQUE.IDENTIFIER", "CATEGORY", "SCIENTIFIC.NAME", "SUBSPECIES.SCIENTIFIC.NAME", "OBSERVATION.COUNT", "COUNTRY.CODE",
      "LATITUDE", "LONGITUDE", "OBSERVATION.DATE", "OBSERVER.ID", "SAMPLING.EVENT.IDENTIFIER", "PROTOCOL.TYPE", "DURATION.MINUTES",
      "EFFORT.DISTANCE.KM", "NUMBER.OBSERVERS", "ALL.SPECIES.REPORTED", "GROUP.IDENTIFIER", "APPROVED")
# Get column indices
inds <- which(names(ebd_sample) %in% vars)
inds

## Test command with smaller dataset
# Desired command
# cut -f1,4,5,6,8,9,14,26,27,28,30,31,32,35,36,38,39,40,42 data/raw/ebd_warblers.csv > data/raw/ebd_warblers_cut.csv
# Create command
bash_command_test <- paste("cut -f", paste(inds, collapse=","), " ../data/ebd/processed/ebd_sample_test.csv > ../data/ebd/processed/rbashtest.csv", sep="")
bash_command_test
# Run command
system(bash_command_test)
# Check result
bash_test <- read_tsv(here("..", "data", "ebd", "processed", "rbashtest.csv"))
head(bash_test)

## Run for full dataset
# Create command
bash_command <- paste("cut -f", paste(inds, collapse=","), " ../data/ebd/processed/ebd_warblers.csv > ../data/ebd/processed/ebd_warblers_cut.csv", sep="")
bash_command
# Run command and calculate runtime
ptm <- proc.time()
system(bash_command)
proc.time() - ptm

#### 5. Test exported file ####
warblers <- read_tsv(here("..", "data", "ebd", "processed", "ebd_warblers_cut.csv"))
names(warblers) = gsub(" ", ".", names(warblers))
head(warblers)
warblers_sampling <- read_tsv(here("..", "data", "ebd", "processed", "ebd_warblers_sampling.csv"))
names(warblers_sampling) = gsub(" ", ".", names(warblers_sampling))
head(warblers_sampling)

# Check summary
warblers_summary <- summary(warblers) # takes some time
warblers_summary

# Check categories
warblers %>%
  filter(CATEGORY == "form") %>%
  select(CATEGORY, SCIENTIFIC.NAME, SUBSPECIES.SCIENTIFIC.NAME) %>%
  head(50)
warblers %>%
  filter(CATEGORY == "intergrade") %>%
  select(CATEGORY, SCIENTIFIC.NAME, SUBSPECIES.SCIENTIFIC.NAME) %>%
  head(50)
warblers %>%
  filter(CATEGORY == "issf") %>%
  select(CATEGORY, SCIENTIFIC.NAME, SUBSPECIES.SCIENTIFIC.NAME) %>%
  head(50)

# Check species
unique(warblers$SCIENTIFIC.NAME) # 63
# Check subspecies
unique(warblers$SUBSPECIES.SCIENTIFIC.NAME) # 41
# Check counts
sort(unique(warblers$OBSERVATION.COUNT)) # X is weird
warblers %>%
  filter(OBSERVATION.COUNT != "X") %>%
  pull(OBSERVATION.COUNT) %>%
  as.numeric %>%
  hist

# Check longitude
filter(warblers, LONGITUDE > 0) # possibly and error

warblers_strict <- warblers %>%
  filter(CATEGORY == "species",
       OBSERVATION.COUNT != "X",
       LONGITUDE < 0,
       PROTOCOL.TYPE == "Traveling",
       APPROVED == 1) %>%
  droplevels
summary(warblers_strict)

# Transform dates
tmp <- head(warblers, 10)
tmp %>%
  mutate_at(vars(OBSERVATION.DATE), list(YEAR = year, MONTH = month, DAY = day))
warblers <- warblers %>%
  mutate_at(vars(OBSERVATION.DATE), list(YEAR = year, MONTH = month, DAY = day))

# Check dates
warblers %>%
  count(YEAR) %>%
  print.data.frame()
warblers %>%
  count(MONTH) %>%
  print.data.frame()
warblers %>%
  count(DAY) %>%
  print.data.frame()

#### 6. Zero-filling ####
# library(auk)
f_ebd <- here("..", "data", "ebd", "processed", "ebd_warblers_cut.csv")
f_sampling <- here("..", "data", "ebd", "processed", "ebd_warblers_sampling.csv")
ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)

## Useful transformations
# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebd_zf <- ebd_zf %>%
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X",
                  NA_character_, observation_count),
                  observation_count = as.integer(observation_count),
                  # effort_distance_km to 0 for non-travelling counts
                  effort_distance_km = if_else(protocol_type != "Traveling",
                  0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started)
  )

## Accounting for variation in detectability (need to change parameters)
# # additional filtering
# ebd_zf_filtered <- ebd_zf %>%
#   filter(
#     # effort filters
#     duration_minutes <= 5 * 60,
#     effort_distance_km <= 5,
#     # last 10 years of data
#     year(observation_date) >= 2009,
#     # 10 or fewer observers
#     number_observers <= 10)

## Select desired variables
ebird <- ebd_zf_filtered %>%
  select(checklist_id, observer_id, sampling_event_identifier,
       scientific_name,
       observation_count, species_observed,
       state_code, locality_id, latitude, longitude,
       protocol_type, all_species_reported,
       observation_date, time_observations_started,
       duration_minutes, effort_distance_km,
       number_observers)

## Save to csv
write_tsv(ebird, here("..", "data", "ebd", "processed", "ebd_warblers_zf.csv"), na = "")
