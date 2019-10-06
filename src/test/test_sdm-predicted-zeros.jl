using Distributed
using JLD2
using ProgressMeter

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
    # Observed coordinates range
    lon_range_obs = extrema(df.longitude)
    lat_range_obs = extrema(df.latitude)
end

## Get the worldclim data
@time wc_vars_pred = pmap(x -> worldclim(x, resolution = "10")[lon_range, lat_range], 1:19);
@time wc_vars_train = pmap(x -> worldclim(x, resolution = "5")[lon_range_obs, lat_range_obs], 1:19);

## Test bioclim & species_bclim functions
occ = warblers_occ[end-4]
pred_vars = wc_vars_pred
train_vars = wc_vars_train

# Predictions
predictions = [bioclim(occ, train_vars[i], train_vars = train_vars[i]) for i in 1:length(pred_vars)];
prediction = reduce(minimum, predictions)
prediction[occ]
# O.0, why???
occ
no_nan = filter(!isnan, prediction[occ])
filter(iszero, prediction[occ])
unique(occ, [:longitude, :latitude])

import SimpleSDMLayers: match_longitude
import SimpleSDMLayers: match_latitude

show(DataFrame(transpose(hcat(map(x -> x[occ[1:5,:]], predictions)...))), allrows=true)
obs = 1
var = 2

observed_values = train_vars[var][occ]
qfinder = ecdf(observed_values)

lon = match_longitude(train_vars[var], occ[obs,:].longitude)
lat = match_latitude(train_vars[var], occ[obs,:].latitude)

train_vars[var].grid[lat, lon] == observed_values[obs]
train_vars[var].grid[lat, lon]
extrema(observed_values)
local_quantile = qfinder(train_vars[var].grid[lat, lon])

grid_lons = map(x -> match_longitude(train_vars[var], x.longitude), eachrow(occ))
grid_lats = map(x -> match_latitude(train_vars[var], x.latitude), eachrow(occ))
grid_pos = DataFrame(grid_lons = grid_lons, grid_lats = grid_lats)
unique(grid_pos, [:grid_lons])
unique(grid_pos, [:grid_lats])
unique(grid_pos, [:grid_lons, :grid_lats])

## Test fixed species_bclim function
test_species = warblers_occ[end-19:end]
pres_abs = map(x -> presence_absence(x, wc_vars_pred[1]), test_species)
@time fixed = @showprogress map(x -> species_bclim(x, wc_vars_pred; train_vars = wc_vars_train), test_species)

filter(!isnan, fixed[end].grid) # species absent
filter(!isnan, fixed[end-1].grid) # species absent
filter(!isnan, fixed[end-2].grid) # species present
filter(!isnan, fixed[1].grid)
npred = map(x -> length(filter(!isnan, x.grid)), fixed)
nobs = map(x -> length(filter(isone, x.grid)), pres_abs)
npred - nobs
(npred - nobs)./nobs .* 100


n = 1
plotSDM(pres_abs[n])
plotSDM(fixed[n])

sp1_pres = presence_absence(warblers_occ[1], wc_vars_pred[1])
sp1_pred = species_bclim(warblers_occ[1], wc_vars_pred; train_vars = wc_vars_train)
fig_sp1_pres = plotSDM(sp1_pres)
fig_sp1_pred = plotSDM(sp1_pred)

savefig(fig_sp1_pred, "test.pdf")
