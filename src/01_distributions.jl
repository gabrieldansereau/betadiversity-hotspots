import Pkg; Pkg.activate(".")
using Distributed
addprocs(9)
@time @everywhere include("src/required.jl")

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "sdm" # desired outcome (required)
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
# save_data = true # should data files be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
  @warn "'outcome' not defined, must be either 'raw' or 'sdm'"
elseif (outcome != "raw" && outcome != "sdm")
  @warn "'outcome' invalid, must be either 'raw' or 'sdm'"
else
  @info "'outcome' currently set to '$(outcome)'"
end

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

## Export distributions
# save_data = true # should data files be overwritten (optional)
if (@isdefined save_data) && save_data == true
    # Export data
    @info "Data exported to file ($(outcome) distributions data)"
    @save "data/jld2/$(outcome)-distributions.jld2" distributions spenames speindex
else
    # Load data
    @info "Data imported from file ($(outcome) distributions data)"
    @load "data/jld2/$(outcome)-distributions.jld2" distributions spenames speindex
end

## Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
sort(pres_counts)

## Plot result
# Species 1
sp1 = "Setophaga_townsendi"
map_sp1 = plotSDM(distributions[speindex[sp1]], c=:BuPu,
                  title = "$(sp1) distribution ($(outcome))",
                  colorbar=:none, dpi=300)
scatter!(map_sp1, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=5)
# Species 2
sp2 = "Setophaga_petechia"
map_sp2 = plotSDM(distributions[speindex[sp2]], c=:BuPu,
                  title = "$(sp2) distribution ($(outcome))",
                  colorbar=:none, dpi=300)
scatter!(map_sp2, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=5)

## Export figures
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) distributions)"
    savefig(map_sp1, "fig/$(outcome)/01_$(outcome)_sp-$(sp1).png")
    savefig(map_sp2, "fig/$(outcome)/01_$(outcome)_sp-$(sp2).png")
else
    @info "Figures not saved ($(outcome) distributions)"
end
