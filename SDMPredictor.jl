struct SDMPredictor{T<:Number}
    grid::Matrix{T}
    left::Float64
    right::Float64
    bottom::Float64
    top::Float64
end

function latitudes(p::SDMPredictor)
    n = size(p.grid, 1)
    grid_size = (p.top-p.bottom)/n/2.0
    centers = range(p.bottom+grid_size; stop=p.top-grid_size, length=n)
    return centers
end

function longitudes(p::SDMPredictor)
    n = size(p.grid, 2)
    grid_size = (p.right-p.left)/n/2.0
    centers = range(p.left+grid_size; stop=p.right-grid_size, length=n)
    return centers
end

function Base.getindex(p::SDMPredictor, longitude::Float64, latitude::Float64)
    i_lon = findmin(abs.(longitude .- longitudes(p)))[2]
    j_lat = findmin(abs.(latitude .- latitudes(p)))[2]
    return p.grid[j_lat, i_lon]
end

function Base.getindex(p::SDMPredictor, longitude::NTuple{2,Float64}, latitude::NTuple{2,Float64})
    m_lon = findmin(abs.(minimum(longitude) .- longitudes(p)))[2]
    M_lon = findmin(abs.(maximum(longitude) .- longitudes(p)))[2]
    m_lat = findmin(abs.(minimum(latitude) .- latitudes(p)))[2]
    M_lat = findmin(abs.(maximum(latitude) .- latitudes(p)))[2]

    n_lon = (p.right-p.left)/size(p.grid, 2)/2.0
    n_left = longitudes(p)[m_lon]-n_lon
    n_right = longitudes(p)[M_lon]+n_lon

    n_lat = (p.top-p.bottom)/size(p.grid, 1)/2.0
    n_bottom = latitudes(p)[m_lat]-n_lat
    n_top = latitudes(p)[M_lat]+n_lat

    n_grid = p.grid[m_lat:M_lat, m_lon:M_lon]
    return SDMPredictor(n_grid, n_left, n_right, n_bottom, n_top)
end

function Base.size(p::SDMPredictor)
    return size(p.grid)
end

function Base.size(p::SDMPredictor, i::Int64)
    return size(p.grid, i)
end
