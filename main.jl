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

# Get some GBIF data
q = Dict{Any,Any}("country" => "CA", "limit" => 100)
occ = occurrences(taxon("Turdus migratorius"), q)
[next!(occ) for i in 1:8]
qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates])

# Get the worldclim data by their layer number
@time wc_vars = [clip(worldclim(i), occ) for i in 1:19];
# Make the prediction for each layer
@time predictions = [bioclim(wc_vars[i], occ) for i in 1:length(wc_vars)];
# Make the final prediction by taking the minimum
@time prediction = reduce(minimum, predictions);

# Filter the predictions so that everything below the 10th percentile is NaN
threshold = quantile(prediction[occ], [0.1])[1]
for i in eachindex(prediction.grid)
    prediction.grid[i] < threshold && (prediction.grid[i] = NaN)
end

# Get the map and plot it
function get_shape(res)
    @assert res ∈ [50,110]
    dir = "https://github.com/nvkelso/natural-earth-vector/raw/master/$(res)m_physical/"
    fn = "ne_$(res)m_land.shp"
    run(`wget $dir/$fn -P /tmp/`)
    handle = open("/tmp/$fn", "r") do io
        read(io, Shapefile.Handle)
    end
    return handle
end
lores = get_shape(50)

sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(1200,600), frame=:box)
xaxis!(sdm_plot, (prediction.left,prediction.right), "Longitude")
yaxis!(sdm_plot, (prediction.bottom,prediction.top), "Latitude")

for p in lores.shapes
    xy = map(x -> (x.x, x.y), p.points)
    sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
    plot!(sdm_plot, sh, c=:lightgrey, lab="")
end

heatmap!(
    sdm_plot,
    longitudes(prediction), latitudes(prediction), prediction.grid, 
    aspectratio=1.3, c=:Oranges, frame=:box
    )

for p in lores.shapes
    xy = map(x -> (x.x, x.y), p.points)
    plot!(sdm_plot, xy, c=:grey, lab="")
end

scatter!(
    sdm_plot,
    longitudes(occ), latitudes(occ),
    c=:black, msw=0.0, ms=0.5, lab=""
    )

savefig("sdm.png")
