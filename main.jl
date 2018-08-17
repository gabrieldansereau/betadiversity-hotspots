using Plots
using Statistics
using GDAL
using GBIF

include("lib/SDMLayer.jl")
include("lib/gdal.jl")
include("lib/worldclim.jl")
include("lib/bioclim.jl")

# Get some GBIF data
q = Dict{Any,Any}("country" => "CA", "limit" => 100)
occ = occurrences(taxon("Lepus townsendii"), q)
[next!(occ) for i in 1:50]
qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates])

# Get the worldclim data by their layer number
@time wc_vars = [clip(worldclim(i), occ) for i in 1:19];
# Make the prediction for each layer
@time predictions = [bioclim(wc_vars[i], occ) for i in 1:length(wc_vars)];
# Make the final prediction by taking the minimum
@time prediction = reduce(minimum, predictions);

heatmap(
        longitudes(prediction), latitudes(prediction), prediction.grid, 
        aspectratio=1.3, c=:Greens
       )
savefig("sdm.png")

scatter!(
         longitudes(occ), latitudes(occ),
         c=:black, msw=0.0, ms=1, lab=""
        )
savefig("sdm_occ.png")
