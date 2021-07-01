## Presence absence distribution functions

"""
    presence_absence(species::AbstractDataFrame, layer::SimpleSDMLayer)

Convert the occurrences the `species` DataFrame to a presence-absence layer
based on the dimensions and coordinates from `layer`.
"""
function presence_absence(species::AbstractDataFrame, layer::SimpleSDMLayer)
    @assert all(in(names(species)), ["longitude", "latitude"]) "Coordinates must be stored in columns called longitude and latitude"

    # Create empty grid for presence-absence data (with NaN)
    distribution_grid = copy(layer.grid)
    replace!(x -> !isnothing(x) ? 0.0 : NaN, distribution_grid)

    # Get unique sites/coordinates
    unique_sites = unique(species, [:longitude, :latitude])

    # Loop for all sites
    lon_all = collect(longitudes(layer))
    lat_all = collect(latitudes(layer))
    for (lon, lat) in zip(unique_sites.longitude, unique_sites.latitude)
        # Get grid position for each site
        i_lon = findmin(abs.(lon .- lon_all))[2]
        j_lat = findmin(abs.(lat .- lat_all))[2]
        # Add 1 per species presence
        distribution_grid[j_lat, i_lon] += 1.0
    end

    # Reduce sites with more than 1 presence to binary value (default)
    replace!(x -> x > 1.0 ? 1.0 : x, distribution_grid)
    # Replace zeros (absences) by nothing
    replace!(distribution_grid, 0.0 => nothing)
    replace!(distribution_grid, NaN => nothing)
    # Create SimpleSDMLayer
    distribution_layer = SimpleSDMResponse(distribution_grid, layer)
    return distribution_layer
end
