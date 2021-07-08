#### Explore eBird Dataset ####
# Mostly following tutorial from http://strimas.com/ebird-best-practices/

source(file.path("src", "required.R"))
library(auk)

# Set ebd path
# auk::auk_set_ebd_path("~/github/data/ebd/raw/")
# Restart session

# Conditional evaluations
# save_prepdata <- TRUE

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
ebird_taxonomy <- as_tibble(ebird_taxonomy)
names(ebird_taxonomy)
ebird_taxonomy

# Extract Parulidae taxonomy information
parulidae_taxonomy <- ebird_taxonomy %>%
  filter(str_detect(family, "Parulidae"))

# Check parulidae taxonomy
distinct(parulidae_taxonomy, family)
parulidae_taxonomy
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
f_ebd <- here("data", "proc", "ebd_sample_test.csv")
auk_filter(ebd_filters, file = f_ebd, overwrite = T)

# Test exported file
test <- read_tsv(f_ebd)
names(test) = gsub(" ", ".", names(test))
head(test)
table(test$SCIENTIFIC.NAME)
table(test$PROTOCOL.TYPE)

#### 3. Extract data from complete database ####

# Check worldwide distribution for a few species (manually from eBird website)
amN <- c("Myioborus miniatus", "Myioborus pictus", "Cardellina rubrifrons")
amS <- c(
  "Myioborus albifrons", "Myioborus melanocephalus",
  "Myioborus ornatus", "Myioborus flavivertex",
  "Myioborus pariae", "Myioborus castaneocapilla",
  "Myioborus brunniceps", "Myiothlypis rivularis",
  "Myiothlypis nigrocristata", "Myiothlypis leucoblephara",
  "Myiothlypis luteoviridis"
  )
amC <- c("Myioborus torquatus", "Cardellina versicolor", "Cardellina rubra")
nul <- c("Myioborus albifacies", "Myioborus cardonai")

# Remove a few species not present in North America (to have 97 at most and run command only once)
to_remove <- c(amS, nul)
fewer_parulidae_species <- parulidae_species[!(parulidae_species %in% to_remove)]
fewer_parulidae_species

# Get ebd file
ebd <- auk_ebd("ebd_relJun-2019.txt", file_sampling = "ebd_sampling_relJun-2019.txt")
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
# (Only if save_prepdata option correctly set)
# save_prepdata <- TRUE
if (exists("save_prepdata") && isTRUE(save_prepdata)) {
  message("Exporting filtered EBD data sets")
  f_ebd <- here("..", "data", "ebd", "processed", "ebd_warblers.csv")
  f_sampling <- here("..", "data", "ebd", "processed", "ebd_warblers_sampling.csv")
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
  auk_filter(ebd_sampling_filters, file = f_sampling)
}
# File is too big, needs to cut in terminal

#### 4. Cut file using bash command ####

# Select variables to keep
vars <- c("GLOBAL.UNIQUE.IDENTIFIER", "CATEGORY", "COMMON.NAME", "SCIENTIFIC.NAME", "SUBSPECIES.SCIENTIFIC.NAME", "OBSERVATION.COUNT", "COUNTRY.CODE",
      "LATITUDE", "LONGITUDE", "OBSERVATION.DATE", "OBSERVER.ID", "SAMPLING.EVENT.IDENTIFIER", "PROTOCOL.TYPE", "DURATION.MINUTES",
      "EFFORT.DISTANCE.KM", "NUMBER.OBSERVERS", "ALL.SPECIES.REPORTED", "GROUP.IDENTIFIER", "APPROVED")

# Get column indices
inds <- which(names(ebd_sample) %in% vars)
inds

# Desired command
# cut -f1,4,5,6,8,9,14,26,27,28,30,31,32,35,36,38,39,40,42 data/raw/ebd_warblers.csv > data/raw/ebd_warblers_cut.csv

# Test command with smaller dataset
bash_command_test <- paste("cut -f", paste(inds, collapse=","), " ./data/proc/ebd_sample_test.csv > ./data/proc/ebd_sample_cut.csv", sep="")
bash_command_test

# Run test command
system(bash_command_test)

# Check result
bash_test <- read_tsv(here("data", "proc", "ebd_sample_cut.csv"))
head(bash_test)

# Create command for full data set
bash_command <- paste("cut -f", paste(inds, collapse=","), " ../data/ebd/processed/ebd_warblers.csv > ../data/ebd/processed/ebd_warblers_cut.csv", sep="")
bash_command

# Run command on full data set and calculate runtime
# (Only if save_prepdata option correctly set)
# save_prepdata <- TRUE
if (exists("save_prepdata") && isTRUE(save_prepdata)) {
  message("Cutting full EBD data set")
  tic("Cutting full EBD data data set")
  system(bash_command)
  toc()
}
