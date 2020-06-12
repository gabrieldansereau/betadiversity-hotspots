## Overloads to existing functions (support for additional types & arguments)

import Base: maximum, minimum, sum, +, -, *, /, min, max
import Statistics: mean, median, std
import SimpleSDMLayers: longitudes
import SimpleSDMLayers: latitudes

## Extract layer values

# Extract layer value from DataFrame (with longitude & latitude columns)
function Base.getindex(p::SimpleSDMLayer, d::DataFrame)
    observations = eltype(p.grid)[]
    for i in 1:nrow(d)
        push!(observations, p[d.longitude[i], d.latitude[i]])
    end
    return observations
end

## Extract all latitudes & longitudes

# Extract longitudes from DataFrame
function longitudes(d::DataFrame)
    l = Float64[]
    for lon in d.longitude
        push!(l, lon)
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

## Overloads for single layers

ops = Symbol.((
    "sum", "maximum", "minimum",
    "mean", "median", "std"
    ))

for op in ops, ty in (:SimpleSDMResponse, :SimpleSDMPredictor)
    eval(quote
        """
            $($op)(l::$($ty){T}) where {T <: Number}
        Applies `$($op)` to an object of type `$($ty)`. This function has been
        automatically generated. Note that this function is only applied to the
        non-`NaN` elements of the layer, and has no method to work on the `dims`
        keyword; the grid itself can be extracted with `convert(Matrix, l)`.
        """
        function $op(l::$ty{T}) where {T <: Number}
            return $op(filter(!isnan, l.grid))
        end
    end)
end

## Overloads for multiple layers

# Math operations
ops_math = Symbol.(("+", "-", "*", "/"))
for op in ops_math
    eval(quote
        function $op(layer1::SimpleSDMLayer, layer2::SimpleSDMLayer)
            SimpleSDMLayers._layers_are_compatible(layer1, layer2)
            newlayer = copy(layer1)
            newlayer.grid = broadcast($op, layer1.grid, layer2.grid)
            return newlayer
        end
    end)
end

# Min/max
for op in (:min, :max)
    eval(quote
        function $op(layer1::SimpleSDMLayer, layer2::SimpleSDMLayer)
            SimpleSDMLayers._layers_are_compatible(layer1, layer2)
            newlayer = copy(layer1)
            for i in eachindex(newlayer)
                newlayer[i] = $op(layer1[i], layer2[i])
            end
            return newlayer
        end
    end)
end

# Mean/std
for op in (:mean, :std)
    eval(quote
        function $op(layers::Array{T}) where {T <: SimpleSDMLayer}
            newlayer = copy(layers[1])
            newlayer.grid = $op(map(x -> x.grid, layers))
            return newlayer
        end
    end)
end
