using Plots
using Statistics
using GDAL
using GBIF

include("lib.jl")

# Extract coordinates
bioclim_codes = lpad.(1:19, 2, "0")
bioclim_variables = [get_tiff_data("assets/wc2.0_bio_10m_$(x).tif") for x in bioclim_codes]

latitudes = collect(range(-90.0; stop=90.0, length=size(bioclim_variables[1], 1)))
longitudes = collect(range(-180.0; stop=180.0, length=size(bioclim_variables[1], 2)))
heatmap(longitudes, latitudes, bioclim_variables[1], aspectratio=1, c=:viridis)

predicates = get_bioclim_values(occ_data, bioclim_variables, longitudes, latitudes)

b_box = [(-155.0, 20.0),(-50.0, 75.0)]

model_per_variable = []
for (i, v) in enumerate(bioclim_variables)
    @info i
    push!(model_per_variable, get_quantile_matrix(v, predicates[:,i], b_box, longitudes, latitudes))
end

retained_models = filter(x -> sum(filter(.!isnan, x)) > 1.90e4, model_per_variable)

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
p1 = heatmap(b_box_lon, b_box_lat, round.(consensus; digits=1), frame=:ticks, c=:viridis, clim=(0,1))
p2 = heatmap(b_box_lon, b_box_lat, round.(consensus; digits=5).>0.0, frame=:ticks, c=:Greys_r, clim=(0,1))
scatter!(p2, obs_lon, obs_lat, leg=false, m=:cross, ms=1, c=:orange, msw=0.0)
for p in (p1, p2)
    xaxis!(p, "Longitude", (minimum(b_box_lon), maximum(b_box_lon)))
    yaxis!(p, "Latitude", (minimum(b_box_lat), maximum(b_box_lat)))
end

p2

(model_per_variable .|> x -> filter(.!isnan, x) |> sum) |> sort |> scatter

using Profile

Profile.clear()
@profile get_quantile_matrix(bioclim_variables[1], predicates[:,1], b_box, longitudes, latitudes)
Juno.profiler()
