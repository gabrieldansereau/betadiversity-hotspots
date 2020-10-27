@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
distributions_old = distributions
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
distributions_new = distributions
distributions = nothing

isequal(distributions_new, distributions_old)
.!map(x -> isequal(distributions_new[x].grid, distributions_old[x].grid), eachindex(distributions_new)) |> findall
inds = .!isequal.(distributions_new[1].grid, distributions_old[1].grid) |> findall
lon_inds = [i[1] for i in inds]
lat_inds = [i[2] for i in inds]
lons = longitudes(distributions_new[1])[lon_inds]
lats = latitudes(distributions_new[1])[lat_inds]

distributions_new[1].grid[inds]
distributions_old[1].grid[inds]
numgrid = fill(NaN, size(distributions_new[1].grid))
numgrid[:] = 1:length(distributions_new[1].grid)
numgrid[inds] # matches the diffed lines from the csv files

### 

function presence_absence_new(species::AbstractDataFrame, copy_layer::SimpleSDMLayer; binary::Bool=true)
    # species: species occurrences in a DataFrame with longitude and latitude columns
    # copy_layer: layer with range of interest and dimensions to copy (not used otherwise)
    # full_range: return species full range, even outside range of interest
    # binary: convert to binary presence-absence values per site

    # Extract coordinates
    coords = (left = copy_layer.left, right = copy_layer.right, bottom = copy_layer.bottom, top = copy_layer.top)
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

function presence_absence_old(species::AbstractDataFrame, copy_layer::SimpleSDMLayer; binary::Bool=true)
    # species: species occurrences in a DataFrame with longitude and latitude columns
    # copy_layer: layer with range of interest and dimensions to copy (not used otherwise)
    # full_range: return species full range, even outside range of interest
    # binary: convert to binary presence-absence values per site

    # Extract coordinates
    coords = (left = copy_layer.left, right = copy_layer.right, bottom = copy_layer.bottom, top = copy_layer.top)
    # Create empty grid for presence-absence data (with NaN)
    distribution_grid = copy(copy_layer.grid)
    replace!(x -> !isnothing(x) ? 0.0 : NaN, distribution_grid)
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

species = warblers[24]
copy_layer = env_vars[1]

testold = presence_absence_old(species, copy_layer)
testnew = presence_absence_new(species, copy_layer)

.!isequal.(testold.grid, testnew.grid) |> findall
# still different

coords = (left = copy_layer.left, right = copy_layer.right, bottom = copy_layer.bottom, top = copy_layer.top)
# Create empty grid for presence-absence data (with NaN)
distribution_grid_old = copy(copy_layer.grid)
replace!(x -> !isnothing(x) ? 0.0 : NaN, distribution_grid_old)
distribution_grid_new = copy(copy_layer.grid)
replace!(x -> !isnothing(x) ? 0.0 : NaN, distribution_grid_new)
# Get unique sites/coordinates
unique_sites = unique(species, [:longitude, :latitude])
# old
for site in eachrow(unique_sites)
    # Get grid position for each site
    i_lon = findmin(abs.(site.longitude .- longitudes(copy_layer)))[2]
    j_lat = findmin(abs.(site.latitude .- latitudes(copy_layer)))[2]
    # Add 1 per species presence
    distribution_grid_old[j_lat, i_lon] += 1.0
end
# new
lon_all = collect(longitudes(copy_layer))
lat_all = collect(latitudes(copy_layer))
for (lon, lat) in zip(unique_sites.longitude, unique_sites.latitude)
    # Get grid position for each site
    i_lon = findmin(abs.(lon .- lon_all))[2]
    j_lat = findmin(abs.(lat .- lat_all))[2]
    # Add 1 per species presence
    distribution_grid_new[j_lat, i_lon] += 1.0
end

isequal(distribution_grid_old, distribution_grid_new)
.!isequal.(distribution_grid_old, distribution_grid_new) |> findall

isequal(lon_all, longitudes(copy_layer))
isequal(lat_all, latitudes(copy_layer))

zipsites = zip(unique_sites.longitude, unique_sites.latitude)
lon, lat = (unique_sites.longitude, unique_sites.latitude)
fill(0., (nrow(unique_sites), 2))
tmp = [[lon, lat] for (lon, lat) in zipsites]
tmp = reduce(hcat, tmp) |> permutedims
isequal(lon, tmp[:,1])
isequal(lat, tmp[:,2])

lon2 = [s.longitude for s in eachrow(unique_sites)]
lat2 = [s.latitude for s in eachrow(unique_sites)]
isequal(lon2, tmp[:,1])
isequal(lat2, tmp[:,2])
isequal(lon, lon2)
isequal(lat, lat2)
isequal(lon, lon2)
isequal(lat, lat2)

isequal.(lon2, tmp[:,1]) |> all
isequal.(lat2, tmp[:,2]) |> all

lon[1] .- lon_all

lon_all == longitudes(copy_layer)
lat_all == latitudes(copy_layer)