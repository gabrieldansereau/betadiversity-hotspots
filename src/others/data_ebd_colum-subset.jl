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

## Speed-up prepare_ebd_data!

# Reset values
cut_file = "./data/raw/ebd_warblers_cut.csv"
cols = [
    "COMMON NAME",
    "SCIENTIFIC NAME",
    "LATITUDE",
    "LONGITUDE",
    "OBSERVATION DATE",
    "GROUP IDENTIFIER",
    "APPROVED",
]

# Load small file
small_df = CSV.read(cut_file, DataFrame; limit=100, select=cols)

# Test replace options (and reload file afterwards)
df = copy(small_df)
@time df."SCIENTIFIC NAME" .= replace.(df."SCIENTIFIC NAME", " " .=> "_"); # 405 alloc
@time df."SCIENTIFIC NAME" = replace.(df."SCIENTIFIC NAME", " " .=> "_"); # 406 alloc
@time [replace(sp, " " => "_") for sp in df."SCIENTIFIC NAME"] # 43.68k alloc
for sp in df."SCIENTIFIC NAME"
    replace!(sp, " " => "_")
end
for d in eachrow(df)
    replace!(d."S", " " => "_")
end
# whatever, leave it as is

## Speed up loading & preparing

# Reset values
cut_file = "./data/raw/ebd_warblers_cut.csv"
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
@time begin
    @time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), DataFrame; select=cols)
    @time prepare_ebd_data!(df)
end
#=
28.041258 seconds (107.50 M allocations: 8.105 GiB, 8.80% gc time, 9.01% compilation time)
67.450549 seconds (112.68 M allocations: 6.743 GiB, 83.60% gc time, 4.39% compilation time)
95.516399 seconds (220.19 M allocations: 14.849 GiB, 61.62% gc time, 5.75% compilation time)
=#

@time begin
    @time small_df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), DataFrame; limit=10, select=cols)
    @time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), DataFrame; select=cols)
    @time prepare_ebd_data!(small_df)
    @time prepare_ebd_data!(df)
end
#=
14.667969 seconds (29.98 M allocations: 1.561 GiB, 3.74% gc time, 0.03% compilation time)
12.707636 seconds (81.45 M allocations: 6.781 GiB, 12.44% gc time, 25.68% compilation time)
 0.669131 seconds (1.01 M allocations: 57.883 MiB, 99.78% compilation time)
38.370923 seconds (112.49 M allocations: 6.733 GiB, 69.06% gc time, 2.13% compilation time)
66.441073 seconds (224.93 M allocations: 15.131 GiB, 43.09% gc time, 7.16% compilation time)
=#

@time begin
    @time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), DataFrame; select=cols)
    @time precompile(prepare_ebd_data!, (DataFrame,))
    @time prepare_ebd_data!(df)
end
#=
28.829118 seconds (107.40 M allocations: 8.099 GiB, 7.73% gc time, 8.81% compilation time)
 0.000006 seconds (4 allocations: 560 bytes)
58.133317 seconds (112.68 M allocations: 6.743 GiB, 80.73% gc time, 1.79% compilation time)
87.018891 seconds (220.10 M allocations: 14.844 GiB, 56.49% gc time, 4.12% compilation time)
=#
# whatever
