## Landcover variables loading functions

struct Copernicus end

import SimpleSDMLayers: SimpleSDMPredictor

"""
    SimpleSDMPredictor(Copernicus, LandCover, layer::Int64; x...)

Return a single landcover layer.
"""
function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Integer=1; resolution::Float64=10.0, path::AbstractString=joinpath("assets", "landcover"), kwargs...)
    ## Get file path
    # List files in path
    lc_files = readdir(path)
    # Filter for selected resolution
    filter!(x -> occursin.("$(Int(resolution))m.tif", x), lc_files)
    # Create path for selected layer only
    p = joinpath.(path, lc_files)[layer]

    ## Load data
    # Read raster dataset
    d = ArchGDAL.read(p)
    # Read raster values from band 1 (only band in this case)
    v = ArchGDAL.read(d, 1)
    # Convert from UInt8 to Float64
    v = convert.(Float64, v)

    ## Preliminary manipulations
    # Reverse rows and permute dimensions
    landcover_mat = permutedims(v[:,end:-1:1])
    # Replace 255 (default no data values) by nothing
    landcover_mat = convert(Array{Union{Float32, Nothing}}, landcover_mat)
    replace!(landcover_mat, 255 => nothing)

    # Fill missing latitudes with nothing (latitude extent is only (-60,80) for landcover data)
    nlat, nlon = size(landcover_mat)
    slim, nlim = abs(-60), 80
    res = Int64(nlat/(nlim+slim))
    south_nothings = fill(nothing, ((90-slim)*res, nlon))
    north_nothings = fill(nothing, ((90-nlim)*res, nlon))
    landcover_grids = vcat(south_nothings, landcover_mat, north_nothings)

    # Convert to SimpleSDMLayers
    landcover_layers = SimpleSDMPredictor(landcover_grids, -180.0, 180.0, -90.0, 90.0)

    return landcover_layers
end

"""
    SimpleSDMPredictor(Copernicus, LandCover, layers::UnitRange{Int64}; x...)

Return a range of landcover layers.
"""
function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layers::AbstractArray; kwargs...)
    @assert eltype(layers) <: Integer
    return [SimpleSDMPredictor(Copernicus, LandCover, l; kwargs...) for l in layers]
end