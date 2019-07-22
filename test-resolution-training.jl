using Plots
using GDAL
using Shapefile
using GBIF
using StatsBase
using Statistics
using DataFrames
using CSV

cd("$(homedir())/github/BioClim/")
include("lib/SDMLayer.jl")
include("lib/gdal.jl")
include("lib/worldclim.jl")
include("lib/bioclim.jl")
include("lib/shapefiles.jl")
include("lib/csvdata.jl")

# Use DataFrame instead of GBIFRecord
df = CSV.read("../data/warblers_qc_2018.csv", header=true, delim="\t")
df = prepare_csvdata(df)
warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]
occ = warblers_occ[1]

# Get the worldclim data by their layer number
@info "Extract and crop bioclim variables"
@time wc_vars_train = [clip(worldclim(i, resolution="5"), occ) for i in 1:19];
@time wc_vars_pred = [clip(worldclim(i, resolution="10"), occ) for i in 1:19];
# Make the prediction for each layer
@info "Predictions for each layer"
@time training = [bioclim_training(wc_vars_train[i], occ) for i in 1:length(wc_vars_train)];
@time predictions = [bioclim_prediction(wc_vars_pred[i], occ, training[i]) for i in 1:length(wc_vars_pred)];

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

sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box,
                title = first(unique(occ.species)))
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

sdm_plot
