using Distributed
using JLD2
using Latexify

@time @everywhere include("src/required.jl")

## EBD data preparation

# Load data from CSV files (from file cut with terminal)
@time df = CSV.read("data/raw/ebd_warblers_cut.csv", header=true, delim="\t")

# Prepare data (arrange values & columns)
df = prepare_ebd_data(df)

# Remove duplicates (BUT not ok for counts, so-so for dates, see group-observation.jl in src/test/)
df_nogroups = filter(x -> ismissing(x[:groupIdentifier]), df)
df_groups = dropmissing(df, :groupIdentifier)
df_groups_unique = unique(df_groups, [:species, :groupIdentifier])
newdf = vcat(df_nogroups, df_groups_unique)

# Select subset with specific columns
select!(newdf, [:species, :year, :latitude, :longitude, :groupIdentifier, :countryCode, :samplingEventIdentifier])
# Remove 1 Aleutian Islands observation with positive longitude
filter!(x -> x.longitude < 0, newdf)

# Clean workspace
df = newdf
df_nogroups = nothing
df_groups = nothing
df_groups_unique = nothing

# Data for table
table1 = by(df, :countryCode, [:countryCode, :species, :samplingEventIdentifier] =>
        x -> (n_obs = length(x.countryCode),
              n_checklist = length(unique(x.samplingEventIdentifier)),
              n_sp = length(unique(x.species))
              ))
tmp = by(df, [:countryCode, :samplingEventIdentifier], nrow)
table2 = by(tmp, :countryCode, [:x1] => x -> (n_sp_moy = mean(x.x1),
                                              n_sp_med = median(x.x1),
                                              n_sp_max = maximum(x.x1)
                                              ))
tables = join(table1, table2, on = :countryCode)

show(md(tables))
latexify(tables) |> print
