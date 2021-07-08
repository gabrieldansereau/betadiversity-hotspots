include("required.jl")

## Conditional arguments
# save_prepdata = true # should prepared data be overwritten (optional)

## EBD data preparation

# Select subset of columns to load from CSV file (overflows memory otherwise)
cols = [
    "COMMON NAME",
    "SCIENTIFIC NAME",
    "LATITUDE",
    "LONGITUDE",
    "OBSERVATION DATE",
    "GROUP IDENTIFIER",
    "APPROVED",
]

# Load data from CSV files (from file cut with terminal)
@time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), DataFrame; select=cols)

# Prepare data (arrange values & columns)
@time prepare_ebd_data!(df)

# Select subset with specific columns
select!(df, [:species, :commonName, :year, :latitude, :longitude, :groupIdentifier])

# Remove 1 Aleutian Islands observation with positive longitude
filter!(:longitude => <(0.0), df)

# Remove duplicates (BUT not ok for counts, so-so for dates, see group-observation.jl in src/test/)
df_nogroups = filter(:groupIdentifier => ismissing, df)
dropmissing!(df, :groupIdentifier)
df_groups_unique = unique(df, [:species, :groupIdentifier])
df = vcat(df_nogroups, df_groups_unique)

# Remove groupIdentifier column
select!(df, Not(:groupIdentifier))

## Export prepared data
# save_prepdata = true # should prepared data be overwritten (optional)
if (@isdefined save_prepdata) && save_prepdata == true
    # Export file
    @info "Data exported to file (data preparation)"
    CSV.write(joinpath("data", "proc", "ebd_warblers_prep.csv"), df; delim="\t")

    # Update placeholder file (as file is too big for version control)
    placeholder_path = joinpath("data", "proc", "ebd_warblers_prep_placeholder.csv")
    open(placeholder_path, "w") do io
        write(io, string(Dates.now()))
    end
else
    @info "Data not exported (data preparation)"
end
