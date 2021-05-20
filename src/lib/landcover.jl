## Landcover variables loading functions

struct Copernicus end

# function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Integer=1; kwargs...)
#     file = _get_raster(EarthEnv, LandCover, layer, full)
#     return geotiff(SimpleSDMPredictor, file; kwargs...)
# end

# function SimpleSDMPredictor(::Type{EarthEnv}, ::Type{LandCover}, layers::AbstractArray; kwargs...)
#     @assert eltype(layers) <: Integer
#     return [SimpleSDMPredictor(EarthEnv, LandCover, l; kwargs...) for l in layers]
# end

import SimpleSDMLayers: SimpleSDMPredictor

# Landcover variables loading
function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layers::Vector{Int64}; resolution::AbstractFloat=10.0, path::AbstractString=joinpath("assets", "landcover/"))
    ## Get file paths
    # List files in path
    lc_files = readdir("$(path)")
    # Filter for selected resolution
    filter!(x -> occursin.("$(Int(resolution))m.tif", x), lc_files)
    # Create paths for selected layers only
    paths = joinpath.(path, lc_files)[layers]

    ## Load data
    # Read raster dataset
    datasets = [ArchGDAL.read(p) for p in paths]
    # Read raster values from band 1 (only band in this case)
    values = [ArchGDAL.read(d, 1) for d in datasets]
    # Convert from UInt8 to Float64
    values_int = [convert.(Float64, v) for v in values]

    ## Preliminary manipulations
    # Reverse rows and permute dimensions
    landcover_mat = [permutedims(v[:,end:-1:1]) for v in values_int]
    # Replace 255 (default no data values) by nothing
    landcover_mat = convert.(Array{Union{Float32, Nothing}}, landcover_mat)
    [replace!(l, 255 => nothing) for l in landcover_mat]

    # Fill missing latitudes with nothing (latitude extent is only (-60,80) for landcover data)
    nlat, nlon = size(landcover_mat[1])
    slim, nlim = abs(-60), 80
    res = Int64(nlat/(nlim+slim))
    south_nothings = fill(nothing, ((90-slim)*res, nlon))
    north_nothings = fill(nothing, ((90-nlim)*res, nlon))
    landcover_grids = [vcat(south_nothings, l, north_nothings) for l in landcover_mat]

    # Convert to SimpleSDMLayers
    landcover_layers = SimpleSDMPredictor.(landcover_grids, -180.0, 180.0, -90.0, 90.0)

    return landcover_layers
end

"""
    SimpleSDMPredictor(Copernicus, LandCover, layer::Int64; x...)

Return a single landcover layer.
"""
SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Int64; x...) = SimpleSDMPredictor(Copernicus, LandCover, [layer]; x...)[1]

"""
    SimpleSDMPredictor(Copernicus, LandCover, layers::UnitRange{Int64}; x...)

Return a range of landcover layers.
"""
SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layers::UnitRange{Int64}; x...) = SimpleSDMPredictor(Copernicus, LandCover, collect(layers); x...)
