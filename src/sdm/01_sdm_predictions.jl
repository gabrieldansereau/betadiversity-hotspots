import Pkg; Pkg.activate(".")
using Distributed
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
# WorldClim data with different training resolutions
@time @everywhere wc_vars_pred = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
@time @everywhere wc_vars_train = pmap(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], [1,12]);
# Landcover data
@time @everywhere lc_vars_pred = landcover(1:10, resolution = "10")[lon_range, lat_range]
@time @everywhere lc_vars_train = landcover(1:10, resolution = "5")[lon_range_obs, lat_range_obs]
# Combine environmental data
@everywhere env_vars_pred = vcat(wc_vars_pred, lc_vars_pred)
@everywhere env_vars_train = vcat(wc_vars_train, lc_vars_train)

## Make predictions for all species
# With different training resolutions
@time predictions = @showprogress pmap(x -> species_bclim(x, pred_vars = env_vars_pred, train_vars = env_vars_train), warblers_occ);

## Export predictions
@save "data/jld2/sdm-predictions-landcover.jld2" predictions

# Test import
@load "data/jld2/sdm-predictions-landcover.jld2" predictions
