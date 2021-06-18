using Pkg: Pkg
Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# Load required results
@load "data/jld2/raw-pres-abs.jld2" pres_abs
@load "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

# Define coordinates range
coords = (left=-145.0, right=-50.0, bottom=20.0, top=75.0)

# Get worldclim variables
@time wc_vars = SimpleSDMPredictor(WorldClim, BioClim, 1:19; resolution=5.0, coords...);

# Get sites pres-abs data
vars_spe = map(x -> vec(x.grid[inds_pred]), pres_abs[1:3])
# Get sites worldclim variables
vars_env = map(x -> vec(x.grid[inds_pred]), wc_vars)

# Get sites latitudes
lats = repeat(collect(latitudes(wc_vars[1])); outer=size(wc_vars[1].grid, 2))[inds_pred]
# Get sites longitudes
lons = repeat(collect(longitudes(wc_vars[1])); inner=size(wc_vars[1].grid, 1))[inds_pred]

# Create dataframe
fluxdf = DataFrame(hcat(lats, lons, vars_spe..., vars_env...), :auto)
# Create column names
fluxnames = [:latitude, :longitude, :sp1, :sp2, :sp3, Symbol.("wcvars", collect(1:19))...]
# Change column names
names!(fluxdf, fluxnames)

# Test coordinates & latitude
test1 = map(x -> wc_vars[1][x.longitude, x.latitude], eachrow(fluxdf))
test2 = fluxdf.wcvars1
test1 .== test2
test1 == test2

# Export dataframe
CSV.write("data/proc/ebd_flux.csv", fluxdf; delim="\t")

# Test load
testdf = CSV.read("data/proc/ebd_flux.csv"; header=true, delim="\t")
