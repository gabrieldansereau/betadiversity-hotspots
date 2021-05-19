struct Copernicus end

function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Integer=1; kwargs...)
    file = _get_raster(EarthEnv, LandCover, layer, full)
    return geotiff(SimpleSDMPredictor, file; kwargs...)
end

function SimpleSDMPredictor(::Type{EarthEnv}, ::Type{LandCover}, layers::AbstractArray; kwargs...)
    @assert eltype(layers) <: Integer
    return [SimpleSDMPredictor(EarthEnv, LandCover, l; kwargs...) for l in layers]
end
