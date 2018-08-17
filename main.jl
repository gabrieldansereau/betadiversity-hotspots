using Plots
using Statistics
using GDAL
using GBIF

# Core type
include("lib/SDMPredictor.jl")

# Read geotiff files
include("lib/gdal.jl")

# Format bioclim data
include("lib/worldclim.jl")

# Other
include("distributions.jl")


predicates = get_bioclim_values(occ_data, bioclim_variables, longitudes, latitudes)

b_box = [(-155.0, 20.0),(-50.0, 75.0)]

model_per_variable = []
for (i, v) in enumerate(bioclim_variables)
    @info i
    push!(model_per_variable, get_quantile_matrix(v, predicates[i], b_box, longitudes, latitudes))
end

retained_models = filter(x -> sum(filter(.!isnan, x)) > 0.0, model_per_variable)

consensus = zeros(Float64, size(retained_models[1]))
for i in eachindex(consensus)
    consensus[i] = minimum(getindex.(retained_models, i))
end

obs_lon = []
obs_lat = []
for r in occ_data
    push!(obs_lon, r.longitude)
    push!(obs_lat, r.latitude)
end

b_box_lon = range(b_box[1][1]; stop=b_box[2][1], length=size(consensus, 2))
b_box_lat = range(b_box[1][2]; stop=b_box[2][2], length=size(consensus, 1))

coordinates = []
for r in occ_data
    push!(coordinates, (r.longitude, r.latitude))
end


raw_suit = [get_value_at_position(coordinate, consensus, collect(b_box_lon), collect(b_box_lat)) for coordinate in coordinates]
suitability = filter(.!isnan, raw_suit)
tr = quantile(suitability, [0.05])[1]

p1 = heatmap(b_box_lon, b_box_lat, round.(consensus; digits=2), frame=:box, c=:viridis, clim=(0,1), leg=false, aspectratio=1.4)

tokeep = zeros(Int64, size(consensus))
for i in eachindex(tokeep)
    if !isnan(consensus[i])
        tokeep[i] = consensus[i] > tr ? 2 : 1
    else
        tokeep[i] = 0
    end
end

p2 = heatmap(b_box_lon, b_box_lat, tokeep, frame=:box, c=:Blues, aspectratio=1.4, leg=false)
scatter!(p2, obs_lon, obs_lat, leg=false, m=:cross, ms=1, c=:orange, msw=0.0)

for p in (p1, p2)
    xaxis!(p, "Longitude", (minimum(b_box_lon), maximum(b_box_lon)))
    yaxis!(p, "Latitude", (minimum(b_box_lat), maximum(b_box_lat)))
end

plot(p1, p2)
