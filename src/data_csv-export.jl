using Distributed
using JLD2

@time @everywhere include("src/required.jl")

## Export data to csv files with names consistent with NEwR

# Load required results
@load "data/jld2/raw-pres-abs.jld2" pres_abs
@load "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

# Define coordinates range
lon_range = (-145.0, -50.0)
lat_range = (20.0, 75.0)

## Extract data
# Get worldclim variables
@time wc_vars = pmap(x -> worldclim(x, resolution = "5")[lon_range, lat_range], 1:19);

# Get presence-absence data for sampled sites, spe matrix
spe = Ypred

# Get environment data for sampled sites
vars_env = map(x -> vec(x.grid[inds_pred]), wc_vars)
# Create env matrix
env = hcat(vars_env...)

# Get sites latitudes
lats = repeat(collect(latitudes(wc_vars[1])), outer=size(wc_vars[1].grid, 2))[inds_pred]
# Get sites longitudes
lons = repeat(collect(longitudes(wc_vars[1])), inner=size(wc_vars[1].grid, 1))[inds_pred]
# Create spa matrix
spa = hcat(lats, lons)

## Convert to dataframes
# Create species dataframe
spedf = DataFrame(spe)
names!(spedf, Symbol.("sp", collect(1:size(spe,2))))
# Create environmental dataframe
envdf = DataFrame(env)
names!(envdf, Symbol.("wc", collect(1:size(env,2))))
# Create species dataframe
spadf = DataFrame(lat = lats, lon = lons)

## Export dataframes
CSV.write("data/proc/pres-abs_spe.csv", spedf, delim="\t")
CSV.write("data/proc/pres-abs_env.csv", envdf, delim="\t")
CSV.write("data/proc/pres-abs_spa.csv", spadf, delim="\t")

# Test load
testdf = CSV.read("data/proc/pres-abs_spe.csv", header=true, delim="\t")
