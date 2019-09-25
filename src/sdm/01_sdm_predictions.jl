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

## Get the worldclim data
@time wc_vars_pred = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);
@time wc_vars_train = pmap(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], 1:19);

## Make predictions for all species
@time predictions = @showprogress pmap(x -> species_bclim(x, wc_vars_pred, train = wc_vars_train), warblers_occ);

## Export predictions
@save "data/jld2/sdm-predictions.jld2" predictions

# Test import
@load "data/jld2/sdm-predictions.jld2" predictions
