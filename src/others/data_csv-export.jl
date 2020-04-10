import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Export data to csv files with names consistent with NEwR

# Load required results
@load "data/jld2/raw-distributions.jld2" distributions spenames speindex
@load "data/jld2/raw-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

# Define coordinates range
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
# Define coordinates range for Quebec (for smaller scale analysis)
coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)

## Get environmental data data
# WorldClim data
@time wc_vars = map(x -> worldclim(x, resolution = "10")[coords], 1:19);
# Landcover data
@time lc_vars = map(x -> landcover(x, resolution = "10")[coords], 1:10);
# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

# Get environment data for sampled sites
env_vars_vec = map(x -> vec(x.grid), env_vars)
# Create env matrix
env_mat = hcat(env_vars_vec...)

# Get sites latitudes
lats = repeat(collect(latitudes(wc_vars[1])), outer=size(wc_vars[1].grid, 2))
# Get sites longitudes
lons = repeat(collect(longitudes(wc_vars[1])), inner=size(wc_vars[1].grid, 1))
# Create spa matrix
spa_mat = hcat(lats, lons)

## Convert to dataframes
# Create species dataframe
spe_full = DataFrame(Y)
rename!(spe_full, Symbol.("sp", 1:ncol(spe_full)))
# Create environmental dataframe
env_full = DataFrame(env_mat)
rename!(env_full, vcat(Symbol.("wc", 1:size(wc_vars, 1)), Symbol.("lc", 1:size(lc_vars, 1))))
# Create species dataframe
spa_full = DataFrame(lat = lats, lon = lons)

## Export dataframes
# CSV.write("data/proc/distributions_spe_full.csv", spe_full, delim="\t") # not needed
CSV.write("data/proc/distributions_env_full.csv", env_full, delim="\t")
CSV.write("data/proc/distributions_spa_full.csv", spa_full, delim="\t")

# Test load
testdf = CSV.read("data/proc/distributions_spa_full.csv", header=true, delim="\t")

## Filter to observed sites ##

# Filter dataframes to observed sites only
spe_obs = spe_full[inds_obs,:]
env_obs = env_full[inds_obs,:]
spa_obs = spa_full[inds_obs,:]

# Export observed dataframes
CSV.write("data/proc/distributions_spe.csv", spe_obs, delim="\t")
CSV.write("data/proc/distributions_env.csv", env_obs, delim="\t")
CSV.write("data/proc/distributions_spa.csv", spa_obs, delim="\t")

#### Restrict data to Quebec ####

## Limit observation range
# Get indexes
inds_qc = findall(x -> (coords_qc.left < x.lon < coords_qc.right) &&
                       (coords_qc.bottom < x.lat < coords_qc.top), eachrow(spa_full))
inds_qcobs = intersect(inds_obs, inds_qc)

# Filter datasets
speqc_full = spe_full[inds_qc,:]
envqc_full = env_full[inds_qc,:]
spaqc_full = spa_full[inds_qc,:]
speqc_obs = spe_full[inds_qcobs,:]
envqc_obs = env_full[inds_qcobs,:]
spaqc_obs = spa_full[inds_qcobs,:]

## Preliminary verifications
# Remove species without observations
cols_obs = findall(x -> sum(x) > 0, eachcol(speqc_obs))
speqc_full = speqc_full[:, cols_obs]
speqc_obs = speqc_obs[:, cols_obs]

## Export QC data
# CSV.write("data/proc/distributions_spe_qc_full.csv", speqc_full, delim="\t") # not needed
CSV.write("data/proc/distributions_env_qc_full.csv", envqc_full, delim="\t")
CSV.write("data/proc/distributions_spa_qc_full.csv", spaqc_full, delim="\t")
CSV.write("data/proc/distributions_spe_qc.csv", speqc_obs, delim="\t")
CSV.write("data/proc/distributions_env_qc.csv", envqc_obs, delim="\t")
CSV.write("data/proc/distributions_spa_qc.csv", spaqc_obs, delim="\t")
