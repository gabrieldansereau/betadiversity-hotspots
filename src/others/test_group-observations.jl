using Distributed
using JLD2

@time @everywhere include("src/required.jl")

df = CSV.read("../data/ebd/processed/ebd_warblers_cut.csv", header=true, delim="\t")

## Data preparation

# Prepare data
df = prepare_ebd_data(df)
# Select specific columns
select!(df, [:species, :groupIdentifier, :observerId, :samplingEventIdentifier, :observationCount, :year, :latitude, :longitude])

# Remove missing groups (unimportant for investigation)
df_gr = dropmissing(df, :groupIdentifier)
# Sort by group & species
df_gp = sort(df_gr, [:groupIdentifier, :species])

#### Investigate differences in group observations

# 1. Groups & observation counts
# 2. Groups & Dates

## 1. Groups & observation counts

# Get unique groups & species combinations
unique_groups = unique(df_gr, [:species, :groupIdentifier])
# Get unique groups, species & counts combinattions
unique_counts = unique(df_gr, [:species, :groupIdentifier, :observationCount])
# Compare number of observations
n_unique_groups = nrow(unique_groups)
n_unique_counts = nrow(unique_counts)
# unique_counts > unique groups, meaning groups don't always have same counts

# Get unique species & groups combination within unique_counts
unique_counts2 = unique(unique_counts, [:species, :groupIdentifier])
n_unique_counts2 = nrow(unique_counts2)
# Compare number of observations
unique_counts2 == unique_groups
# equal, confirms it's only groups that don't have the same count

# Get non unique observations
nonu = unique_counts[nonunique(unique_counts, [:species, :groupIdentifier]),:]
sort!(nonu, [:groupIdentifier, :species])
# Doesn't work!!!
# nonunique() returns rows equal to one of the previous rows, hence not the first duplicate

# Alternative way: using groupby for similar observations
grouped = groupby(unique_counts, [:groupIdentifier, :species])
# Get non unique observations as groups with > 1 observation
nonu_grouped = grouped[nrow.(grouped) .> 1]

nonu_grouped[1]
nonu_grouped[2]
nonu_grouped[3]
nonu_grouped[4]
nonu_grouped[100]
# Confirms observation counts are NOT the same per observer
# Also, observerId can come twice within a group obs.
# I guess one is the master observation and the 2nd second was edited to be the observer's personal count
# Hence, right thing to do seems to be: remove duplicates BUT select observation with highest count

## 2. Groups & dates

# Get unique groups, species & dates combinations
unique_dates = unique(df_gr, [:species, :groupIdentifier, :observationDate])
n_unique_dates = nrow(n_unique_dates)
# Compare number of observations
n_unique_dates - n_unique_groups

# Get non unique observations (groups-species-dates)
grouped_dates1 = groupby(unique_dates, [:groupIdentifier, :species])
nonu_grouped_dates1 = grouped_dates1[nrow.(grouped_dates1) .> 1]
# Reformat to dataframe & select columns
dates_test1 = vcat(nonu_grouped_dates1...)[[:species, :groupIdentifier, :observerId, :observationDate, :observationCount]]
# Check results
nrow(dates_test1) # 7895 observations, 3940 groups
show(first(dates_test1,100), allrows=true)
# Most differences seem to be minor or typos (± a few days), but sometimes the difference is a few years
# If the year is ever important for the analyses, those should be removed

# Get non unique dates from dataset of unique observation counts (less important, but was my initial intent)
grouped_dates = map(x -> unique(x, :observationDate), nonu_grouped)
nonu_grouped_dates = grouped_dates[nrow.(grouped_dates) .> 1]
nonu_grouped_dates[1] |> DataFrame |> x -> select(x, [:species, :groupIdentifier, :observerId, :observationDate])
nonu_grouped_dates[2] |> DataFrame |> x -> select(x, [:species, :groupIdentifier, :observerId, :observationDate])
nonu_grouped_dates[3] |> DataFrame |> x -> select(x, [:species, :groupIdentifier, :observerId, :observationDate])
nonu_grouped_dates[4] |> DataFrame |> x -> select(x, [:species, :groupIdentifier, :observerId, :observationDate])
test = vcat(nonu_grouped_dates...)[[:species, :groupIdentifier, :observerId, :observationDate]]
nrow(test) # 146 observations, 73 groups
show(test, allrows=true)
