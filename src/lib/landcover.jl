## Landcover variables loading functions

struct Copernicus end

import SimpleSDMLayers: SimpleSDMPredictor

"""
    SimpleSDMPredictor(Copernicus, LandCover, layer::Int64; x...)

Return a single landcover layer.
"""
function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Integer=1; resolution::Float64=10.0, path::AbstractString=joinpath("assets", "landcover"), kwargs...)
    @assert resolution in [5.0, 10.0]
    1 ≤ layer ≤ 10 || throw(ArgumentError("The layer must be between 1 and 10"))

    # List files in path
    lc_files = readdir(path)
    # Filter for selected resolution
    filter!(x -> occursin.("$(Int(resolution))m.tif", x), lc_files)
    # Create path for selected layer only
    p = joinpath.(path, lc_files)[layer]
    
    # Read layer
    layer = geotiff(SimpleSDMPredictor, p; kwargs...)
    # Convert from UInt8 to Float32
    layer = convert(Float32, layer)

    return layer
end

"""
    SimpleSDMPredictor(Copernicus, LandCover, layers::UnitRange{Int64}; x...)

Return a range of landcover layers.
"""
function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layers::AbstractArray; kwargs...)
    @assert eltype(layers) <: Integer
    return [SimpleSDMPredictor(Copernicus, LandCover, l; kwargs...) for l in layers]
end