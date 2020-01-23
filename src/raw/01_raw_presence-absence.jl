import Pkg; Pkg.activate(".")
using Distributed
addprocs(4)
@time @everywhere include("src/required.jl")

outcome = "raw"

## Get & prepare data
@time begin
    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)

    # Load data from CSV files
    df = CSV.read("data/proc/ebd_warblers_prep.csv", header=true, delim="\t")
    df = filter(x -> lon_range[1] < x[:longitude] < lon_range[2], df)
    filter!(x -> lat_range[1] < x[:latitude] < lat_range[2], df)

    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

end

## Get environmental data data
# WorldClim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
# Landcover data
@time lc_vars = landcover(1:10, resolution = "10")[lon_range, lat_range]
# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

## Create function to convert occurrence to presence-absence based on a SimpleSDMLayer
@everywhere function presence_absence(species::DataFrame, copy_layer::SimpleSDMLayer; binary::Bool=true)
    # Create empty grid for presence-absence data (with NaN)
    distributions_grid = copy(copy_layer.grid)
    replace!(x -> !isnan(x) ? 0.0 : x, distributions_grid)
    # Get unique sites/coordinates
    unique_sites = unique(species, [:longitude, :latitude])
    # Loop for all sites
    for site in eachrow(unique_sites)
        # Get grid position
        i_lon = findmin(abs.(site.longitude .- longitudes(copy_layer)))[2]
        j_lat = findmin(abs.(site.latitude .- latitudes(copy_layer)))[2]
        # Add 1 per species presence
        distributions_grid[j_lat, i_lon] += 1.0
    end
    # Check sites with presence
    filter(x -> x > 0.0, distributions_grid)
    # Reduce sites with more than 1 presence
    if binary == true
        replace!(x -> x > 1.0 ? 1.0 : x, distributions_grid)
    end
    # Create SimpleSDMLayer
    distributions_layer = SimpleSDMResponse(distributions_grid,
                              copy_layer.left, copy_layer.right,
                              copy_layer.bottom, copy_layer.top)
    return distributions_layer
end
# Loop function for each species
@time distributions = @showprogress pmap(x -> presence_absence(x, env_vars[1]), warblers_occ)
# @time distributions1 = pmap(x -> presence_absence(x, env_vars[1], binary=false), warblers_occ)
# Export result
@save "data/jld2/$(outcome)-pres-abs.jld2" distributions
@load "data/jld2/$(outcome)-pres-abs.jld2" distributions

# Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
sort(pres_counts)

## Plot result
map_sp1 = plotSDM(distributions[1], c=:BuPu)
title!(map_sp1, "$(warblers_occ[1].species[1]) distribution (presence-absence)")
map_sp2 = plotSDM(distributions[13], c=:BuPu)
title!(map_sp2, "$(warblers_occ[13].species[1]) distribution (presence-absence)")

# Export figure
#=
savefig(map_sp1, "fig/$(outcome)/01_$(outcome)_sp-$(warblers_occ[1].species[1]).pdf")
savefig(map_sp2, "fig/$(outcome)/01_$(outcome)_sp-$(warblers_occ[13].species[1]).pdf")
=#
