import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

## Get & prepare data
# Load data from CSV files
@time df = CSV.read(joinpath("data", "proc", "ebd_warblers_prep.csv"), header=true, delim="\t")
# Separate species
warblers = [df[df.species .== u,:] for u in unique(df.species)]
# Reorder species by frequency
sort!(warblers, by = x -> nrow(x), rev = true)
# Extract species names
spenames = [w.species[1] for w in warblers]
# Create index Dict for species names
speindex = indexmap(spenames)

# Define coordinates range
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
# Observed coordinates range
coords_obs = (left = minimum(df.longitude), right = maximum(df.longitude),
              bottom = minimum(df.latitude), top = maximum(df.latitude))

## Get environmental data (with different training resolutions)
# WorldClim data
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, [1, 12]; resolution = 10.0, coords...);
# Landcover data
lc_vars = map(x -> SimpleSDMPredictor(Copernicus, LandCover, x; resolution = 10.0)[coords], 1:10);
# Training data with finer resolution
wc_vars_train = SimpleSDMPredictor(WorldClim, BioClim, [1, 12]; resolution = 5.0, coords_obs...);
lc_vars_train = map(x -> SimpleSDMPredictor(Copernicus, LandCover, x; resolution = 5.0)[coords_obs], 1:10);

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)
env_vars_train = vcat(wc_vars_train, lc_vars_train)

## Partition training-validation datasets
import MLJ.partition
# Create empty sets
train_sets = Array{DataFrame, 1}()
valid_sets = Array{DataFrame, 1}()
for w in warblers
    # Get training and validation indices
    train_idx, valid_idx = partition(eachindex(eachrow(w)), 0.7, shuffle = true, rng = 42)
    # Extract training and validations sets
    push!(train_sets, w[train_idx, :])
    push!(valid_sets, w[valid_idx, :])
end
train_sets
valid_sets

## Get distribution for all species
# Get sdm distributions (with different training resolutions)
@time distributions = @showprogress pmap(x -> bioclim(x, env_vars, training_layers = env_vars_train), train_sets[[24,3]]);

## Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
sort(pres_counts)

## Plot result
# Species 1
sp1 = "Setophaga_townsendi"
map_sp1 = plotSDM(distributions[1], c=:BuPu,
                  title = "$(sp1) distribution",
                  colorbar=:none, dpi=300)
scatter!(map_sp1, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=5)
# Species 2
sp2 = "Setophaga_petechia"
map_sp2 = plotSDM(distributions[speindex[sp2]], c=:BuPu,
                  title = "$(sp2) distribution",
                  colorbar=:none, dpi=300)
scatter!(map_sp2, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=5)
