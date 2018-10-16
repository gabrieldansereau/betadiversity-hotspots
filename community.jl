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

function gbifdata(sp)
    @info sp
    q = Dict{Any,Any}("limit" => 200, "country" => "CA")
    occ = occurrences(sp, q)
    [next!(occ) for i in 1:29]
    qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates])
    return occ
end

parulidae = taxon("Parulidae"; strict=false)
parulidae_latest_200 = occurrences(parulidae, Dict("limit"=>200, "country"=>"CA"))
canadian_warblers = unique([p.taxon for p in parulidae_latest_200])

warblers_occ = gbifdata.(canadian_warblers)
lon_range = (-136.0, -58.0)
lat_range = (40.5, 56.0)

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
