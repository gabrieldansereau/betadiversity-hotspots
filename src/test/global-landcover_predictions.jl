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
oldlc_vars_train = map(x -> landcover(x, resolution = "5")[lon_range, lat_range], 1:10);

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)
env_vars_train = vcat(wc_vars_train, lc_vars_train)
oldenv_vars_train = vcat(wc_vars_train, oldlc_vars_train)

## BIOCLIM
occ = warblers[36]
newres = bioclim(occ, env_vars, train_vars = env_vars_train)
oldres = bioclim(occ, env_vars, train_vars = oldenv_vars_train)
notrainres = bioclim(occ, env_vars)

plotSDM(newres) # empty collection
plotSDM(oldres)
plotSDM(notrainres)

# Apply 1st part of BIOCLIM on each environmental variable
newpredictions = [bioclim_singlevar(occ, env_vars[i], train_vars = env_vars_train[i]) for i in 1:length(env_vars)];
oldpredictions = [bioclim_singlevar(occ, env_vars[i], train_vars = oldenv_vars_train[i]) for i in 1:length(env_vars)];
# Reduce to single layer with minimum values
newprediction = reduce(minimum, newpredictions);
oldprediction = reduce(minimum, oldpredictions);
# Replace zeros by NaN
replace!(x -> x == 0.0 ? NaN : x, newprediction.grid);
replace!(x -> x == 0.0 ? NaN : x, oldprediction.grid);

filter(x -> x > 0, newprediction.grid)
filter(x -> x > 0, oldprediction.grid)

map(l -> length(filter(x -> x > 0, l.grid)), newpredictions) # problem with vars 6-8
map(l -> length(filter(x -> x > 0, l.grid)), oldpredictions)

filter(x -> x > 0, bioclim_singlevar(occ, env_vars[6], train_vars = env_vars_train[6]).grid)

## BIOCLIM singlevar
# occ: occurences of a single species as a DataFrame with latitude and longitude columns
# pred_vars: environmental variables used for prediction
# train_vars: optional, training environmental variables, can use a different resolution

# Get observed environmental values (training values)
newtrain_vars = env_vars_train[6]
oldtrain_vars = oldenv_vars_train[6]
newobserved_values = newtrain_vars[occ]
oldobserved_values = oldtrain_vars[occ]

filter(x -> x > 0, newobserved_values)
filter(x -> x > 0, oldobserved_values)
filter(!iszero, newobserved_values)
filter(!iszero, oldobserved_values)

# Create ECDF function to extract quantile value
newqfinder = ecdf(newobserved_values)
oldqfinder = ecdf(oldobserved_values)
# Create empty array for local quantile values
newlq = zeros(Float64, size(env_vars[6]))
oldlq = zeros(Float64, size(env_vars[6]))
# Loop for all sites (prediction values)
for i in eachindex(env_vars[6].grid)
    if isnan(env_vars[6].grid[i])
        # Set value to NaN if env value is already NaN
        newlocal_quantile = NaN
        oldlocal_quantile = NaN
    else
        # Get quantile rank if not NaN
        newlocal_quantile = newqfinder(env_vars[6].grid[i])
        oldlocal_quantile = oldqfinder(env_vars[6].grid[i])
        # Replace values greater than 0.5 to set both tails equal
        if newlocal_quantile > 0.5
            newlocal_quantile = 1.0-newlocal_quantile
        end
        if oldlocal_quantile > 0.5
            oldlocal_quantile = 1.0-oldlocal_quantile
        end
        # Scale back between 0 and 1
        newlocal_quantile = 2.0 * newlocal_quantile
        oldlocal_quantile = 2.0 * oldlocal_quantile
    end
    # Collect quantile values
    newlq[i] = newlocal_quantile
    oldlq[i] = oldlocal_quantile
end
filter(x -> x > 0, newlq)
filter(x -> x > 0, oldlq)

##

# Small internal function for BIOCLIM. Ensures values are equal for both tails and between 0-1
bcscore(x) = isnan(x) ? NaN : (x > 0.5 ? 2*(1-x) : 2x)
# BIOCLIM prediction for single variable
function bioclimpred(layer::T, df::D; train_layer::T = layer) where {T <: SimpleSDMLayer, D <: DataFrame}
  obs = train_layer[df]
  filter!(!isnan, obs)
  pred = similar(layer)
  if length(unique(obs)) <= 1
      pred.grid = replace(x -> x == obs[1] ? 1.0 : 0.0, layer.grid)
  else
      qf = ecdf(obs)
      pred.grid = bcscore.(qf.(layer.grid))
  end
  for idx in findall(isnan, layer.grid)
    pred.grid[idx] = NaN
  end
  return pred
end

newtest = bioclimpred(env_vars[4], occ, train_layer = env_vars_train[4])
oldtest = bioclimpred(env_vars[6], occ, train_layer = oldenv_vars_train[6])

unique(newtest.grid)
unique(oldtest.grid)

layer = env_vars[6]
df = occ
train_layer = oldenv_vars_train[6]
