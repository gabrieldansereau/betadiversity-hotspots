using Distributed
using JLD2

# required.jl equivalent
begin
    using Plots
    using GDAL
    using GBIF
    using StatsBase
    using Statistics
    using DataFrames
    using CSV
    using Dates
end

# Add SimpleSDMLayers.jl
#=
using Pkg
Pkg.develop(PackageSpec(path="$(homedir())/github/SimpleSDMLayers.jl"))
=#
using SimpleSDMLayers

# Change directory
cd("../BioClim/")

## Get & prepare data
@time begin
    # Load data from CSV files
    df = CSV.read("data/proc/ebd_warblers_prep.csv", header=true, delim="\t")
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
end

## worldclim.jl
# worldclim() --> default resolution is different
test1 = pmap(worldclim, 1:2)
test1b = worldclim(1:2) # new way for multiple layers
solu1 = pmap(x -> worldclim(x, resolution = "5"), 1:2)
size(test1[1])
size(test1b[1])
size(solu1[1])

test1 = solu1
temp = test1[1]

## shapefiles.jl
# missing from SimpleSDMLayer???

#### SDMLayer.jl
## basics.jl
getindex(temp, 1, 1)
temp[1,1] # doesn't work
temp[1:1, 1:1] # returns as layer
temp[1:1, 1:1].grid[1] # shitty way of doing it
temp.grid[1,1] # works for both

getindex(temp, -180.0, -90.0)

## Latitude-longitude inversion problem
getindex(temp, -89.0, -180.0) # NaN, outside of grid
getindex(temp, -180.0, -89.0) # correct value
temp[-90.0, -180.0] # NaN, outside of grid
temp[-180.0, -89.0] # correct value
temp.grid[1, 12] # incorrect value
temp.grid[12, 1] # correct value, inverted call
temp[(-180.0, -178.0), (-90.0, -89.0)].grid # 12 latitude x 24 longitude array
temp.grid[1:24, 1:12] # incorrect, 24 latitude x 12 longitude array
temp.grid[1:12, 1:24] # correct, inverted call

## Get out of package development mode
]free ~/github/SimpleSDMLayer.jl
