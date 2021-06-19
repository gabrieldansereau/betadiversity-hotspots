## Presence absence distribution functions

## Convert occurrences to presence-absence layer based on an existing SimpleSDMLayer
function presence_absence(
    species::AbstractDataFrame, copy_layer::SimpleSDMLayer; binary::Bool=true
)
    # species: species occurrences in a DataFrame with longitude and latitude columns
    # copy_layer: layer with range of interest and dimensions to copy (not used otherwise)
    # full_range: return species full range, even outside range of interest
    # binary: convert to binary presence-absence values per site

    # Create empty grid for presence-absence data (with NaN)
    distribution_grid = copy(copy_layer.grid)
    replace!(x -> !isnothing(x) ? 0.0 : NaN, distribution_grid)
    # Get unique sites/coordinates
    unique_sites = unique(species, [:longitude, :latitude])
    # Loop for all sites
    lon_all = collect(longitudes(copy_layer))
    lat_all = collect(latitudes(copy_layer))
    for (lon, lat) in zip(unique_sites.longitude, unique_sites.latitude)
        # Get grid position for each site
        i_lon = findmin(abs.(lon .- lon_all))[2]
        j_lat = findmin(abs.(lat .- lat_all))[2]
        # Add 1 per species presence
        distribution_grid[j_lat, i_lon] += 1.0
    end
    # Reduce sites with more than 1 presence to binary value (default)
    if binary == true
        replace!(x -> x > 1.0 ? 1.0 : x, distribution_grid)
    end
    # Replace zeros (absences) by nothing
    replace!(distribution_grid, 0.0 => nothing)
    replace!(distribution_grid, NaN => nothing)
    # Create SimpleSDMLayer
    distribution_layer = SimpleSDMResponse(distribution_grid, copy_layer)
    return distribution_layer
end
