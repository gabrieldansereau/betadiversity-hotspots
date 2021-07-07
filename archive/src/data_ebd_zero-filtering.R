#### EBD Zero-filling ####
# Mostly following tutorial from http://strimas.com/ebird-best-practices/
# !!! CAREFUL !!! Overflows memory on full data set

source(file.path("src", "required.R"))
library(auk)
library(lubridate)

## Prepare data
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
