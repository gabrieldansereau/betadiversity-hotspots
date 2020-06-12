import Base: maximum, minimum, sum, +, -, *, /, min, max
import Statistics: mean, median, std

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
