using Distributed
using JLD2
addprocs(9)

@time @everywhere include("src/required.jl")

## Get & prepare data
@time begin
    # Load data from CSV files
    df = CSV.read("../data/ebd/ebd_warblers_cut.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_ebd_data(df)
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x)[lon_range, lat_range], 1:19);

## Make predictions for all species
@time predictions = pmap(x -> species_bclim(x, wc_vars), warblers_occ);

## Export predictions
@save "../data/predictions-ebd.jld2" predictions

# Test import
@load "../data/predictions-ebd.jld2" predictions
