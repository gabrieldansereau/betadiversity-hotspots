## Presence absence distribution functions

## Convert occurrences to presence-absence layer based on an existing SimpleSDMLayer
function presence_absence(species::DataFrame, copy_layer::SimpleSDMLayer; full_range::Bool=false, binary::Bool=true)
    # species: species occurrences in a DataFrame with longitude and latitude columns
    # copy_layer: layer with range of interest and dimensions to copy (not used otherwise)
    # full_range: return species full range, even outside range of interest
    # binary: convert to binary presence-absence values per site

    # Extract coordinates
    coords = (left = copy_layer.left, right = copy_layer.right, bottom = copy_layer.bottom, top = copy_layer.top)
    # Filter observations to range of interest (default)
    if full_range == false
        filter!(x -> coords.left < x[:longitude] < coords.right, species)
        filter!(x -> coords.bottom < x[:latitude] < coords.top, species)
    end
    # Create empty grid for presence-absence data (with NaN)
    distribution_grid = copy(copy_layer.grid)
    replace!(x -> !isnan(x) ? 0.0 : x, distribution_grid)
    # Get unique sites/coordinates
    unique_sites = unique(species, [:longitude, :latitude])
    # Loop for all sites
    for site in eachrow(unique_sites)
        # Get grid position for each site
        i_lon = findmin(abs.(site.longitude .- longitudes(copy_layer)))[2]
        j_lat = findmin(abs.(site.latitude .- latitudes(copy_layer)))[2]
        # Add 1 per species presence
        distribution_grid[j_lat, i_lon] += 1.0
    end
    # Replace zeros (absences) by NaN
    replace!(distribution_grid, 0.0 => NaN)
    # Reduce sites with more than 1 presence to binary value (default)
    if binary == true
        replace!(x -> x > 1.0 ? 1.0 : x, distribution_grid)
    end
    # Create SimpleSDMLayer
    distribution_layer = SimpleSDMResponse(distribution_grid,
                              copy_layer.left, copy_layer.right,
                              copy_layer.bottom, copy_layer.top)
    return distribution_layer
end
