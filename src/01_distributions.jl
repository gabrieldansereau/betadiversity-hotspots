import Pkg; Pkg.activate(".")
using Distributed
addprocs(9)
@time @everywhere include("src/required.jl")

# outcome = "sdm"

## Get & prepare data
# Load data from CSV files
@time df = CSV.read("data/proc/ebd_warblers_prep.csv", header=true, delim="\t")
# Separate species
warblers = [df[df.species .== u,:] for u in unique(df.species)]
# Reorder species by frequency
sort!(warblers, by = x -> nrow(x), rev = true)
# Extract species names
spenames = [w.species[1] for w in warblers]
# Create index Dict for species names
speindex = indexmap(spenames)

# Define coordinates range
lon_range = (-145.0, -50.0)
lat_range = (20.0, 75.0)
# Observed coordinates range
lon_range_obs = extrema(df.longitude)
lat_range_obs = extrema(df.latitude)

## Get environmental data (with different training resolutions)
# WorldClim data
@time wc_vars = map(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
# Landcover data
@time lc_vars = pmap(x -> landcover(x, resolution = "10")[lon_range, lat_range], 1:10);
# Training data with finer resolution
if outcome == "raw"
    # Set resolution to 10 # CAN'T BE FINER FOR RAW ANALYSES FOR NOW
    @time wc_vars_train = map(x -> worldclim(x, resolution = "10")[lon_range_obs, lat_range_obs], [1,12]);
    @time lc_vars_train = map(x -> landcover(x, resolution = "10")[lon_range, lat_range], 1:10)
elseif outcome == "sdm"
    # Set resolution to 5
    @time wc_vars_train = map(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], [1,12]);
    @time lc_vars_train = map(x -> landcover(x, resolution = "5")[lon_range, lat_range], 1:10)
end

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)
env_vars_train = vcat(wc_vars_train, lc_vars_train)

## Get distribution for all species
if outcome == "raw"
    # Get raw distributions
    @time distributions = @showprogress pmap(x -> presence_absence(x, env_vars[1]), warblers)
    # @time distributions = @showprogress pmap(x -> presence_absence(x, env_vars_train[1], full_range = true, binary = false), warblers)
elseif outcome == "sdm"
    # Get sdm distributions (with different training resolutions)
    @time distributions = @showprogress pmap(x -> species_bclim(x, env_vars, train_vars = env_vars_train), warblers);
end

## Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
sort(pres_counts)

## Export distributions
@save "data/jld2/$(outcome)-distributions.jld2" distributions spenames speindex
# Test import
@load "data/jld2/$(outcome)-distributions.jld2" distributions spenames speindex

## Plot result
sp1 = "Setophaga_townsendi"
map_sp1 = plotSDM(distributions[speindex[sp1]], c=:BuPu)
title!(map_sp1, "$(sp1) townsendi distribution ($(outcome))")
sp2 = "Setophaga_petechia"
map_sp2 = plotSDM(distributions[speindex[sp2]], c=:BuPu)
title!(map_sp2, "$(sp2) distribution ($(outcome))")

# Export figure
#=
savefig(map_sp1, "fig/$(outcome)/01_$(outcome)_sp-$(sp1).pdf")
savefig(map_sp2, "fig/$(outcome)/01_$(outcome)_sp-$(sp2).pdf")
=#
