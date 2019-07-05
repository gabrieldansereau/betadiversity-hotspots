using Plots
using GDAL
using Shapefile
using GBIF
using StatsBase
using Statistics
using JLD2
using FileIO
using Dates

cd("$(homedir())/github/BioClim/")
include("lib/SDMLayer.jl")
include("lib/gdal.jl")
include("lib/worldclim.jl")
include("lib/bioclim.jl")
include("lib/shapefiles.jl")

function gbifdata(sp)
    @info sp
    q = Dict{Any,Any}("limit" => 200, "country" => "CA")
    occ = occurrences(sp, q)
    [next!(occ) for i in 1:29]
    qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates])
    return occ
end

parulidae = taxon("Parulidae"; strict=false)
q = Dict("limit"=>200, "country"=>"CA")
q[String(:family)*"Key"] = getfield(parulidae, :family).second
parulidae_latest_200 = occurrences(q)
canadian_warblers = unique([p.taxon for p in parulidae_latest_200])

canadian_warblers = canadian_warblers[.!ismissing.(getfield.(canadian_warblers, :species))]

# @elapsed warblers_occ = gbifdata.(canadian_warblers)
# @save "../data/warblers_gbifdata.jld2" warblers_occ
@load "../data/warblers_gbifdata.jld2" warblers_occ

lon_range = (-136.0, -58.0)
lat_range = (40.5, 56.0)

# using CSV
# warblers_occ = CSV.read("$(homedir())/github/data/warblers_mtl.csv", header=true, delim="\t")

@time wc_vars = [worldclim(i)[lon_range, lat_range] for i in 1:19];

function species_bclim(occ, vars)
    predictions = [bioclim(vars[i], occ) for i in 1:length(vars)];
    prediction = reduce(minimum, predictions);
    for i in eachindex(prediction.grid)
        if prediction.grid[i] == 0.0
            prediction.grid[i] = NaN
        end
    end
    return prediction
end

@time predictions = [species_bclim(w, wc_vars) for w in warblers_occ]

# Essai ecdf
bioclim(wc_vars[1], warblers_occ[1])

observed_values = wc_vars[1][warblers_occ[1]]
qfinder = ecdf(observed_values)
local_quantile = [qfinder(wc_vars[1].grid[i]) for i in eachindex(wc_vars[1].grid)]
lq = zeros(Float64, size(wc_vars[1]))
for i in eachindex(wc_vars[1].grid)
    if isnan(wc_vars[1].grid[i])
        local_quantile = NaN
    else
        local_quantile = qfinder(wc_vars[1].grid[i])
        if local_quantile > 0.5
            local_quantile = 1.0-local_quantile
        end
        local_quantile = 2.0 * local_quantile
    end
    lq[i] = local_quantile
end
lq
prediction = SDMLayer(lq, wc_vars[1].left, wc_vars[1].right, wc_vars[1].bottom, wc_vars[1].top)
heatmap(prediction.grid)

# Try to plot each variable (for 1 species at a time)
plot_array = Any[]
for i in 1:9
    prediction = bioclim(wc_vars[i], warblers_occ[1])
    push!(plot_array, heatmap(prediction.grid))
end
plot(plot_array...)

function pielou(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ sum(A)
    return abs(sum(p.*log.(p))/length(p))
end

function shannon(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ sum(A)
    return abs(sum(p.*log.(p)))
end

output = zeros(Float64, size(predictions[1]))
@time for i in 1:size(output, 1), j in 1:size(output, 2)
    x = getindex.(predictions, i, j)
    output[i,j] = shannon(x)
end

evenness = SDMLayer(output, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)
test = SDMLayer([[1 2]; [3 4]], -160.0, 160.0, -80.0, 80.0)

worldmap = clip(worldshape(50), evenness)

sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box)
xaxis!(sdm_plot, (evenness.left,evenness.right), "Longitude")
yaxis!(sdm_plot, (evenness.bottom,evenness.top), "Latitude")

for p in worldmap
    sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
    plot!(sdm_plot, sh, c=:lightgrey, lab="")
end

heatmap!(
    sdm_plot,
    longitudes(evenness), latitudes(evenness), evenness.grid,
    aspectratio=92.60/60.75, c=:BuPu,
    clim=(0.0, maximum(filter(!isnan, evenness.grid)))
    )

for p in worldmap
    xy = map(x -> (x.x, x.y), p.points)
    plot!(sdm_plot, xy, c=:grey, lab="")
end

savefig("warblers.png")
