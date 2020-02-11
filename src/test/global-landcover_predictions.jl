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

plotSDM(newres)
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
