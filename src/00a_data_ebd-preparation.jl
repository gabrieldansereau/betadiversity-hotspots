import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# save_prepdata = true # should prepared data be overwritten (optional)

## EBD data preparation

# Load data from CSV files (from file cut with terminal)
@time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), header=true, delim="\t")

# Prepare data (arrange values & columns)
df = prepare_ebd_data(df)

# Remove duplicates (BUT not ok for counts, so-so for dates, see group-observation.jl in src/test/)
df_nogroups = filter(x -> ismissing(x[:groupIdentifier]), df)
df_groups = dropmissing(df, :groupIdentifier)
df_groups_unique = unique(df_groups, [:species, :groupIdentifier])
newdf = vcat(df_nogroups, df_groups_unique)

# Select subset with specific columns
select!(newdf, [:species, :year, :latitude, :longitude, :groupIdentifier])
# Remove 1 Aleutian Islands observation with positive longitude
filter!(x -> x.longitude < 0, newdf)

## Export prepared data
# save_prepdata = true # should prepared data be overwritten (optional)
if (@isdefined save_prepdata) && save_prepdata == true
    @info "Data exported to file (data preparation)"
    CSV.write(joinpath("data", "proc", "ebd_warblers_prep.csv"), newdf, delim="\t")
else
    @info "Data not exported (data preparation)"
end
