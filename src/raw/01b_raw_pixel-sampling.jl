Pkg.activate(".")
using Distributed
addprocs(4)
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

## Create function to find site position in layer grid
function siteinlayer(site, layer)
    # Get grid position
    i_lon = findmin(abs.(site.longitude .- longitudes(layer)))[2]
    j_lat = findmin(abs.(site.latitude .- latitudes(layer)))[2]
    return (j_lat, i_lon)
end
## Create function to convert occurrence to presence-absence based on a SimpleSDMLayer
@everywhere function presence_absence(species::DataFrame, copy_layer::SimpleSDMLayer; binary::Bool=true)
# Create empty grid for presence-absence data (with NaN)
pres_abs_grid = copy(copy_layer.grid)
replace!(x -> !isnan(x) ? 0.0 : x, pres_abs_grid)
# Get unique sites/coordinates
#=
unique_sites = unique(species, [:longitude, :latitude])
unique_sites.i_lon = 0
unique_sites.j_lat = 0
# Loop for all sites
for site in eachrow(unique_sites)
# Get grid position
site.i_lon = findmin(abs.(site.longitude .- longitudes(copy_layer)))[2]
site.j_lat = findmin(abs.(site.latitude .- latitudes(copy_layer)))[2]
end
unique(unique_sites, [:i_lon, :j_lat])
=#
@time grid_pos = [siteinlayer(obs, wc_vars[1]) for obs in eachrow(df)]
@time grid_pos = pmap(x -> siteinlayer(x, wc_vars[1]), eachrow(df))

# Check sites with presence
filter(x -> x > 0.0, pres_abs_grid)
# Reduce sites with more than 1 presence
if binary == true
    replace!(x -> x > 1.0 ? 1.0 : x, pres_abs_grid)
end
# Create SimpleSDMLayer
pres_abs_layer = SimpleSDMResponse(pres_abs_grid,
                          copy_layer.left, copy_layer.right,
                          copy_layer.bottom, copy_layer.top)
    return pres_abs_layer
end
# Loop function for each species
@time pres_abs = pmap(x -> presence_absence(x, wc_vars[1]), warblers_occ)
# @time pres_abs1 = pmap(x -> presence_absence(x, wc_vars[1], binary=false), warblers_occ)

# Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in pres_abs]
sort(pres_counts) # one with 0 sites

## Plot result
plotSDM(pres_abs[1])
