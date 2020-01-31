import SimpleSDMLayers: longitudes
import SimpleSDMLayers: latitudes

function Base.getindex(p::SimpleSDMLayer, r::GBIFRecord)
    return p[r.longitude, r.latitude]
end

function Base.getindex(p::SimpleSDMLayer, r::GBIFRecords)
    observations = eltype(p.grid)[]
    for record in r
        push!(observations, p[record])
    end
    return observations
end
function Base.getindex(p::SimpleSDMLayer, d::DataFrame)
    observations = eltype(p.grid)[]
    for i in 1:length(d.species)
        push!(observations, p[d.longitude[i], d.latitude[i]])
    end
    return observations
end
function clip(p::SimpleSDMLayer, r::Union{GBIFRecords,DataFrame})
    lats = latitudes(r)
    lons = longitudes(r)
    return p[(minimum(lons)-1.0, maximum(lons)+1.0), (minimum(lats)-1.0, maximum(lats)+1.0)]
end


function longitudes(r::GBIFRecords)
    l = Float64[]
    for record in r
        push!(l, record.longitude)
    end
    return l
end
function longitudes(d::DataFrame)
    l = Float64[]
    for lon in d.longitude
        push!(l, lon)
    end
    return l
end

function latitudes(r::GBIFRecords)
    l = Float64[]
    for record in r
        push!(l, record.latitude)
    end
    return l
end
function latitudes(d::DataFrame)
    l = Float64[]
    for lat in d.latitude
        push!(l, lat)
    end
    return l
end

function Base.minimum(p1::SimpleSDMLayer, p2::SimpleSDMLayer)
    n1 = copy(p1.grid)
    for i in eachindex(p1.grid)
        n1[i] = min(p1.grid[i], p2.grid[i])
    end
    return SimpleSDMResponse(n1, p1.left, p1.right, p1.bottom, p1.top)
end
