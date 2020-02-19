import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Export data to csv files with names consistent with NEwR

# Load required results
@load "data/jld2/raw-distributions.jld2" distributions spenames speindex
@load "data/jld2/raw-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

# Define coordinates range
lon_range = (-145.0, -50.0)
lat_range = (20.0, 75.0)

## Get environmental data data
# WorldClim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);
# Landcover data
@time lc_vars = map(x -> landcover(x, resolution = "10")[lon_range, lat_range], 1:10);
# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

# Get presence-absence data for sampled sites, spe matrix
spe = Yobs

# Get environment data for sampled sites
vars_env = map(x -> vec(x.grid[inds_obs]), env_vars)
# Create env matrix
env = hcat(vars_env...)

# Get sites latitudes
lats = repeat(collect(latitudes(env_vars[1])), outer=size(env_vars[1].grid, 2))[inds_obs]
# Get sites longitudes
lons = repeat(collect(longitudes(env_vars[1])), inner=size(env_vars[1].grid, 1))[inds_obs]
# Create spa matrix
spa = hcat(lats, lons)

## Convert to dataframes
# Create species dataframe
spedf = DataFrame(spe)
rename!(spedf, Symbol.("sp", collect(1:size(spe,2))))
# Create environmental dataframe
envdf = DataFrame(env)
rename!(envdf, vcat(Symbol.("wc", collect(1:size(wc_vars,1 ))), Symbol.("lc", collect(1:size(lc_vars, 1)))))
# Create species dataframe
spadf = DataFrame(lat = lats, lon = lons)

## Export dataframes
CSV.write("data/proc/distributions_spe.csv", spedf, delim="\t")
CSV.write("data/proc/distributions_env.csv", envdf, delim="\t")
CSV.write("data/proc/distributions_spa.csv", spadf, delim="\t")

# Test load
testdf = CSV.read("data/proc/distributions_spe.csv", header=true, delim="\t")
