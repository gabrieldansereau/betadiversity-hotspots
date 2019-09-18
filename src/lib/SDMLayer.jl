import SimpleSDMLayers: longitudes
import SimpleSDMLayers: latitudes

#=
struct SimpleSDMLayer{T<:Number}
    grid::Matrix{T}
    left::Float64
    right::Float64
    bottom::Float64
    top::Float64
end

function Base.size(p::SimpleSDMLayer)
    return size(p.grid)
end

function Base.size(p::SimpleSDMLayer, i::Int64)
    return size(p.grid, i)
end

function Base.stride(p::SimpleSDMLayer)
    lat_stride = (p.top-p.bottom)/size(p, 1)/2.0
    lon_stride = (p.right-p.left)/size(p, 2)/2.0
    return (lon_stride, lat_stride)
end

function Base.stride(p::SimpleSDMLayer, i::Int64)
    @assert i ∈ [1,2]
    return stride(p)[i]
end

function latitudes(p::SimpleSDMLayer)
    grid_size = stride(p, 2)
    centers = range(p.bottom+grid_size; stop=p.top-grid_size, length=size(p, 1))
    return centers
end

function longitudes(p::SimpleSDMLayer)
    grid_size = stride(p, 1)
    centers = range(p.left+grid_size; stop=p.right-grid_size, length=size(p, 2))
    return centers
end

function Base.eachindex(p::SimpleSDMLayer)
    return eachindex(p.grid)
end

function Base.getindex(p::SimpleSDMLayer, i::Int64, j::Int64)
    return p.grid[i,j]
end

function Base.getindex(p::SimpleSDMLayer, i::Int64)
    return p.grid[i]
end

function Base.getindex(p::SimpleSDMLayer, longitude::Float64, latitude::Float64)
    longitude < p.left && return NaN
    longitude > p.right && return NaN
    latitude < p.bottom && return NaN
    latitude > p.top && return NaN
    i_lon = findmin(abs.(longitude .- longitudes(p)))[2]
    j_lat = findmin(abs.(latitude .- latitudes(p)))[2]
    return p[j_lat, i_lon]
end

function Base.getindex(p::SimpleSDMLayer, i::UnitRange{Int64}, j::UnitRange{Int64})
    return p.grid[i,j]
end

function Base.getindex(p::SimpleSDMLayer, longitude::NTuple{2,Float64}, latitude::NTuple{2,Float64})
    m_lon = findmin(abs.(minimum(longitude) .- longitudes(p)))[2]
    M_lon = findmin(abs.(maximum(longitude) .- longitudes(p)))[2]
    m_lat = findmin(abs.(minimum(latitude) .- latitudes(p)))[2]
    M_lat = findmin(abs.(maximum(latitude) .- latitudes(p)))[2]

    n_lon = (p.right-p.left)/size(p, 2)/2.0
    n_left = longitudes(p)[m_lon]-n_lon
    n_right = longitudes(p)[M_lon]+n_lon

    n_lat = (p.top-p.bottom)/size(p, 1)/2.0
    n_bottom = latitudes(p)[m_lat]-n_lat
    n_top = latitudes(p)[M_lat]+n_lat

    n_grid = p[m_lat:M_lat, m_lon:M_lon]
    return SimpleSDMLayer(n_grid, n_left, n_right, n_bottom, n_top)
end
=#

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
    return SimpleSDMLayer(n1, p1.left, p1.right, p1.bottom, p1.top)
end
