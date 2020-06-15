## Overloads to existing functions (support for additional types & arguments)

import Base: getindex, maximum, minimum, sum, +, -, *, /, min, max, broadcast, broadcast!
import Statistics: mean, median, std
import SimpleSDMLayers: longitudes, latitudes

## DataFrame overloads

# Extract layer values from DataFrame (with longitude & latitude columns)
function getindex(layer::SimpleSDMLayer, d::DataFrame)
    observations = eltype(layer.grid)[]
    for i in 1:nrow(d)
        push!(observations, layer[d.longitude[i], d.latitude[i]])
    end
    return observations
end

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

# Clip layer to extent of DataFrame occurrences
function clip(layer::SimpleSDMLayer, r::DataFrame)
    lats = latitudes(r)
    lons = longitudes(r)
    clip_coords = (left = minimum(lons)-1.0, right = maximum(lons)+1.0,
                   bottom = minimum(lats)-1.0, top = maximum(lats)+1.0)
    return layer[clip_coords]
end

## Overloads for single layers

# Basic operations
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

# Extend broadcast on grid values, returning layer
function broadcast(f, layer::SimpleSDMLayer, As...)
    newgrid = broadcast(f, layer.grid, As...)
    newlayer = SimpleSDMResponse(newgrid, layer.left, layer.right, layer.bottom, layer.top)
    return newlayer
end

## Overloads for multiple layers

# Math operations
ops_math = Symbol.(("+", "-", "*", "/"))
for op in ops_math
    eval(quote
        function $op(layer1::SimpleSDMLayer, layer2::SimpleSDMLayer)
            SimpleSDMLayers._layers_are_compatible(layer1, layer2)
            newgrid = broadcast($op, layer1.grid, layer2.grid)
            newlayer = SimpleSDMResponse(newgrid, layer1.left, layer1.right, layer1.bottom, layer1.top)
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
            newlayer = convert(SimpleSDMResponse, newlayer)
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
            newgrid = $op(map(x -> x.grid, layers))
            newlayer = SimpleSDMResponse(newgrid, layers[1].left, layers[1].right, layers[1].bottom, layers[1].top)
            return newlayer
        end
    end)
end
