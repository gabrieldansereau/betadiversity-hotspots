#### Select column subset from the EBD ####
# Because ebd_warblers_cut.csv cannot be loaded in memory anymore

include("../required.jl")

# Set path to files
cut_file = "./data/raw/ebd_warblers_cut.csv"
head_file = "./data/proc/ebd_warblers_cut_head.csv"
new_file = "./data/proc/ebd_warblers_recut.csv"

# Create head csv files to get column names
# cut_cmd = `head $cut_file \> $head_file`
# run(cut_cmd) # pipelines don't work directly in Julia
run(pipeline(`head $cut_file`, head_file))

# Test read
@time head1 = CSV.read(head_file, DataFrame)

# Even simpler: load first rows only
@time head2 = CSV.read(cut_file, DataFrame; limit=9)
isequal(head1, head2)

# Select subset of names
show(stdout, "text/plain", names(head2))
names_subset = [
    # "GLOBAL UNIQUE IDENTIFIER",
    # "CATEGORY",
    "COMMON NAME",
    "SCIENTIFIC NAME",
    # "SUBSPECIES SCIENTIFIC NAME",
    # "OBSERVATION COUNT", # needed in tests
    # "COUNTRY CODE",
    "LATITUDE",
    "LONGITUDE",
    "OBSERVATION DATE",
    # "OBSERVER ID", # needed in tests
    # "SAMPLING EVENT IDENTIFIER", # needed in tests
    # "PROTOCOL TYPE",
    # "DURATION MINUTES",
    # "EFFORT DISTANCE KM",
    # "NUMBER OBSERVERS",
    # "ALL SPECIES REPORTED",
    "GROUP IDENTIFIER",
    "APPROVED",
]
findall(in(names_subset), names(head2))

# Load names subset on head lines first (to test)
head3 = CSV.read(cut_file, DataFrame; limit=9, select=names_subset)

# Test prepare_ebd_data! function on selected columns
prepare_ebd_data!(head3)
select!(head3, [:species, :commonName, :year, :latitude, :longitude, :groupIdentifier])

# Load names subset on full data set (to see if it fit in memory)
@time df = CSV.read(cut_file, DataFrame; select=names_subset)
# Fits in memory 🎉

# Test prepare_ebd_data! function on selected columns
prepare_ebd_data!(df)
select!(df, [:species, :commonName, :year, :latitude, :longitude, :groupIdentifier])
# Works 🎉
