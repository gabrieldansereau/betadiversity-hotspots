using Distributed
using JLD2
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
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "5")[lon_range, lat_range], 1:19);

## Make predictions for all species
@time predictions = @showprogress pmap(x -> species_bclim(x, wc_vars), warblers_occ);

## Export predictions
@save "data/jld2/predictions-ebd.jld2" predictions

# Test import
@load "data/jld2/predictions-ebd.jld2" predictions
