using Plots
using GDAL
using Shapefile
using GBIF
using StatsBase
using Statistics

include("lib/SDMLayer.jl")
include("lib/gdal.jl")
include("lib/worldclim.jl")
include("lib/bioclim.jl")
include("lib/shapefiles.jl")

# Get some GBIF data
q_orig = Dict{Any,Any}("limit" => 200)
occ_orig = occurrences(taxon("Cypripedium reginae"), q)
q = Dict{Any,Any}("limit" => 200, "scientificName" => "Cypripedium reginae", "hasCoordinate" => true)
occ= occurrences(q)
[next!(occ) for i in 1:19]
function is_ca_or_us(r::GBIFRecord)
    r.countryCode ∈ ["CA", "US"]
end
qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates, is_ca_or_us])

# Get the worldclim data by their layer number
@info "Extract and crop bioclim variables"
@time wc_vars = [clip(worldclim(i), occ) for i in 1:19];
# Make the prediction for each layer
@info "Predictions for each layer"
@time predictions = [bioclim(wc_vars[i], occ) for i in 1:length(wc_vars)];
# Make the final prediction by taking the minimum
@info "Minimum-consensus aggregation"
@time prediction = reduce(minimum, predictions);
# Get the threshold for NaN given a percentile
@info "Threshold estimation"
@time threshold = first(quantile(prediction[occ], [0.05]))
@info "5% threshold:\t$(round(threshold; digits=3))"
# Filter the predictions based on the threshold
@info "Final prediction filtering"
@time for i in eachindex(prediction.grid)
    prediction.grid[i] < threshold && (prediction.grid[i] = NaN)
end

worldmap = clip(worldshape(50), prediction)

sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box)
xaxis!(sdm_plot, (prediction.left,prediction.right), "Longitude")
yaxis!(sdm_plot, (prediction.bottom,prediction.top), "Latitude")

for p in worldmap
    sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
    plot!(sdm_plot, sh, c=:lightgrey, lab="")
end

heatmap!(
    sdm_plot,
    longitudes(prediction), latitudes(prediction), prediction.grid,
    aspectratio=92.60/60.75, c=:BuPu,
    clim=(0.0, maximum(filter(!isnan, prediction.grid)))
    )

for p in worldmap
    xy = map(x -> (x.x, x.y), p.points)
    plot!(sdm_plot, xy, c=:grey, lab="")
end

#=
scatter!(
    sdm_plot,
    longitudes(occ), latitudes(occ),
    c=:black, msw=0.0, ms=0.1, lab="",
    alpha=0.5
    )
=#

savefig("sdm.png")
