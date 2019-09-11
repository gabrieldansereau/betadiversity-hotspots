using Distributed
using JLD2
addprocs(4)

@time @everywhere include("src/required.jl")

## Get & prepare data
@time begin
    # Load data from CSV files
    df = CSV.read("../data/ebd/ebd_warblers_prep.csv", header=true, delim="\t")
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x)[lon_range, lat_range], 1:19);

## Create function to convert occurrence to presence-absence based on a SDMLayer
@everywhere function presence_absence(species::DataFrame, copy_layer::SDMLayer; binary::Bool=true)
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
    # Create SDMLayer
    pres_abs_layer = SDMLayer(pres_abs_grid,
                              copy_layer.left, copy_layer.right,
                              copy_layer.bottom, copy_layer.top)
    return pres_abs_layer
end
# Loop function for each species
@time pres_abs = pmap(x -> presence_absence(x, wc_vars[1]), warblers_occ)
# @time pres_abs1 = pmap(x -> presence_absence(x, wc_vars[1], binary=false), warblers_occ)
# Export result
@save "../data/pres-abs-ebd.jld2" pres_abs
@load "../data/pres-abs-ebd.jld2" pres_abs

# Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in pres_abs]
sort(pres_counts) # one with 0 sites

## Plot result
plotSDM(pres_abs[1])

# Get dimensions
nsites = prod(size(pres_abs[1]))
nspecies = length(pres_abs)
# Create Y
Y = zeros(Int64, (nsites, nspecies))
# Fill Y with community predictions
@progress for gc in 1:nsites # loop for all sites
    # Group predictions for all species in site
    R = map(x -> x.grid[gc], pres_abs)
    # Fill Y with binary values
    global Y[gc,:] = isone.(R)
end
# Export matrix Y
@save "../data/Y-pres-abs-ebd.jld2" Y
@load "../data/Y-pres-abs-ebd.jld2" Y

## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
## Option 2: Calculate LCBD only for sites with predictions
# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
# Select sites with predictions only
Ypred = Y[inds_pred,:]
# Compute BD statistics
resBDpred = BD(Ypred)
# Extract LCBD values
LCBDi = resBDpred.LCBDi
# Scale LCBDi values to maximum value
LCBDi = LCBDi./maximum(LCBDi)

## Arrange LCBD values as grid
# Create empty grid
t_lcbd = fill(NaN, size(pres_abs[1]))
# Fill in grid
t_lcbd[inds_pred] = LCBDi
# Create SDMLayer with LCBD values
LCBD = SDMLayer(t_lcbd, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot results
lcbd_plot = plotSDM(LCBD, type="lcbd")
title!(lcbd_plot, "LCBD values per site (relative to maximum)")

## Save result
# savefig(lcbd_plot, "fig/pres-abs-ebd.pdf")
