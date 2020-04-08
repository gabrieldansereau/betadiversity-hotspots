import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

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
wc_vars = map(x -> worldclim(x, resolution = "10")[lon_range, lat_range], [1,12]);
# Landcover data
lc_vars = map(x -> landcover(x, resolution = "10")[lon_range, lat_range], 1:10);
# Training data with finer resolution
wc_vars_train = map(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], [1,12]);
lc_vars_train = map(x -> landcover(x, resolution = "5")[lon_range_obs, lat_range_obs], 1:10);

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)
env_vars_train = vcat(wc_vars_train, lc_vars_train)

## Get distribution for all species
# Runs only if create_distributions = true, as it can take a while. Loaded from file otherwise
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
if (@isdefined create_distributions) && create_distributions == true
    @info "$(outcome) distributions to be created"
    # Select function to run given desired outcome
    if outcome == "raw"
        # Get raw distributions
        @time distributions = @showprogress pmap(x -> presence_absence(x, env_vars[1]), warblers)
        # @time distributions = @showprogress pmap(x -> presence_absence(x, env_vars_train[1], full_range = true, binary = false), warblers)
    elseif outcome == "sdm"
        # Get sdm distributions (with different training resolutions)
        @time distributions = @showprogress pmap(x -> bioclim(x, env_vars, training_layers = env_vars_train), warblers);
    end
end
