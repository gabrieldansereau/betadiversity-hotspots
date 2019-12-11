using Distributed
using JLD2
using ProgressMeter
addprocs(9)

@time @everywhere include("src/required.jl")

## Get & prepare data
@time begin
    # Load data from CSV files
    df = CSV.read("data/proc/ebd_warblers_prep.csv", header=true, delim="\t")
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
    # Observed coordinates range
    lon_range_obs = extrema(df.longitude)
    lat_range_obs = extrema(df.latitude)
end

## Get environmental data
# WorldClim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
# WorldClim data with different training resolutions
#=
@time wc_vars_pred = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);
@time wc_vars_train = pmap(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], 1:19);
=#
# Landcover data
@time lc_vars = load_landcover(lon_range, lat_range)
# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

## Make predictions for all species
# With environmental data
@time predictions = @showprogress pmap(x -> species_bclim(x, env_vars), warblers_occ);
# With different training resolutions
#=
@time predictions = @showprogress pmap(x -> species_bclim(x, wc_vars_pred, train_vars = wc_vars_train), warblers_occ);
=#

## Export predictions
@save "data/jld2/sdm-predictions-landcover.jld2" predictions

# Test import
@load "data/jld2/sdm-predictions-landcover.jld2" predictions
