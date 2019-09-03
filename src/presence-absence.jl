using Distributed
using JLD2

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

## Create new layer
pres_abs = SDMLayer{Float64}[]
@time for species in warblers_occ
    # Create empty grid for presence-absence data (with NaN)
    pres_abs_grid = copy(wc_vars[1].grid)
    replace!(x -> !isnan(x) ? 0.0 : x, pres_abs_grid)
    # Reduce to unique sites/coordinates
    unique_sites = unique(species, [:longitude, :latitude])
    # Loop for all sites
    for site in eachrow(unique_sites)
        # Get grid position
        i_lon = findmin(abs.(site.longitude .- longitudes(wc_vars[1])))[2]
        j_lat = findmin(abs.(site.latitude .- latitudes(wc_vars[1])))[2]
        # Add 1 per species presence
        pres_abs_grid[j_lat, i_lon] += 1.0
    end
    # Check sites with presence
    filter(!iszero, pres_abs_grid)
    # Reduce sites with more than 1 presence
    replace!(x -> x > 1.0 ? 1.0 : x, pres_abs_grid)
    # Create SDMLayer
    pres_abs_layer = SDMLayer(pres_abs_grid,
                                wc_vars[1].left, wc_vars[1].right,
                                wc_vars[1].bottom, wc_vars[1].top)
    # Export result
    push!(pres_abs, pres_abs_layer)
end

# Get dimensions
nsites = prod(size(predictions[1]))
nspecies = length(predictions)
# Create Y
Y = zeros(Int64, (nsites, nspecies))
# Fill Y with community predictions
@progress for gc in eachindex(predictions[1].grid) # loop for all sites
    # Group predictions for all species in site
    R = map(x -> x.grid[gc], predictions)
    # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
    global Y[gc,:] = .!isnan.(R)
end
