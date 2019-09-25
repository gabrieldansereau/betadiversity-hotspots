using Distributed
using JLD2
addprocs(4)

@time @everywhere include("src/required.jl")

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

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x, resolution = "5")[lon_range, lat_range], 1:19);

## Create function to convert occurrence to presence-absence based on a SimpleSDMLayer
@everywhere function presence_absence(species::DataFrame, copy_layer::SimpleSDMLayer; binary::Bool=true)
    # Create empty grid for presence-absence data (with NaN)
    pres_abs_grid = copy(copy_layer.grid)
    replace!(x -> !isnan(x) ? 0.0 : x, pres_abs_grid)
    # Get unique sites/coordinates
    unique_sites = unique(species, [:longitude, :latitude])
    # Loop for all sites
    for site in eachrow(unique_sites)
        # Get grid position
        i_lon = findmin(abs.(site.longitude .- longitudes(copy_layer)))[2]
        j_lat = findmin(abs.(site.latitude .- latitudes(copy_layer)))[2]
        # Add 1 per species presence
        pres_abs_grid[j_lat, i_lon] += 1.0
    end
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
using ProgressMeter
@time pres_abs = @showprogress pmap(x -> presence_absence(x, wc_vars[1]), warblers_occ)
# @time pres_abs1 = pmap(x -> presence_absence(x, wc_vars[1], binary=false), warblers_occ)
# Export result
@save "data/jld2/raw-pres-abs.jld2" pres_abs
@load "data/jld2/raw-pres-abs.jld2" pres_abs

# Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in pres_abs]
sort(pres_counts)

## Plot result
map_sp1 = plotSDM(pres_abs[1])
title!(map_sp1, "$(warblers_occ[1].species[1]) distribution (presence-absence)")
map_sp2 = plotSDM(pres_abs[13])
title!(map_sp2, "$(warblers_occ[13].species[1]) distribution (presence-absence)")

# Export figure
#=
savefig(map_sp1, "fig/raw/01_raw_sp-$(warblers_occ[1].species[1]).pdf")
savefig(map_sp2, "fig/raw/01_raw_sp-$(warblers_occ[13].species[1]).pdf")
=#
