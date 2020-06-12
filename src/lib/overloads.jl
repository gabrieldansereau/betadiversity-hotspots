## Overloads to existing functions (support for additional types & arguments)

import SimpleSDMLayers: longitudes
import SimpleSDMLayers: latitudes

## Extract layer values

# Extract layer value from single GBIFRecord
function Base.getindex(p::SimpleSDMLayer, r::GBIFRecord)
    return p[r.longitude, r.latitude]
end

# Extract layer value from multiple GBIFRecords
function Base.getindex(p::SimpleSDMLayer, r::GBIFRecords)
    observations = eltype(p.grid)[]
    for record in r
        push!(observations, p[record])
    end
    return observations
end

# Extract layer value from DataFrame (with longitude & latitude columns)
function Base.getindex(p::SimpleSDMLayer, d::DataFrame)
    observations = eltype(p.grid)[]
    for i in 1:nrow(d)
        push!(observations, p[d.longitude[i], d.latitude[i]])
    end
    return observations
end

## Extract all latitudes & longitudes

# Extract longitudes from GBIFRecords
function longitudes(r::GBIFRecords)
    l = Float64[]
    for record in r
        push!(l, record.longitude)
    end
    return l
end

# Extract longitudes from DataFrame
function longitudes(d::DataFrame)
    l = Float64[]
    for lon in d.longitude
        push!(l, lon)
    end
    return l
end

# Extract latitudes from GBIFRecords
function latitudes(r::GBIFRecords)
    l = Float64[]
    for record in r
        push!(l, record.latitude)
    end
    return l
end

# Extract latitudes from DataFrame
function latitudes(d::DataFrame)
    l = Float64[]
    for lat in d.latitude
        push!(l, lat)
    end
    return l
end

## Other layer manipulations

# Clip layer to DataFrame/GBIFRecords occurrences extent
function clip(p::SimpleSDMLayer, r::Union{GBIFRecords,DataFrame})
    lats = latitudes(r)
    lons = longitudes(r)
    return p[(minimum(lons)-1.0, maximum(lons)+1.0), (minimum(lats)-1.0, maximum(lats)+1.0)]
end
