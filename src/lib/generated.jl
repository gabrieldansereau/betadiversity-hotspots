import Base: maximum, minimum, sum, +, -, *, /
import Statistics: mean, median, std

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

# Test layers
layer1 = prob_distrib[1]
layer2 = prob_distrib[2]

# Test addition
+(layer1.grid, layer2.grid)
+(layer1, layer2)
+(layer1, layer2).grid
layer1 + worldclim(1)

# Test on array of layers
prob_distrib[1] + prob_distrib[2] + prob_distrib[3]
tmp = reduce(+, prob_distrib)
plotSDM2(tmp, c = :viridis)
histogram(tmp)
sort(unique(tmp.grid))

plotSDM2(reduce(+, lower_distrib), c = :viridis)
plotSDM2(reduce(+, upper_distrib), c = :viridis)

uncertainty = upper_distrib .- lower_distrib
plotSDM2(reduce(+, uncertainty))
plotSDM2(reduce(mean, uncertainty)) # not working

grids = map(x -> x.grid, uncertainty)
heatmap(mean(grids)) # yeah ok, need to implement elementwise mean

# Test multiplication
layer1.grid * layer2.grid # matrix product, not working
layer1.grid .* layer2.grid # elementwise product
.*(layer1.grid, layer2.grid)
broadcast(*, layer1.grid, layer2.grid)
@which broadcast(*, layer1.grid, layer2.grid)
plotSDM2(reduce(*, prob_distrib)) # elementwise layer product

# Test division
layer1.grid / layer2.grid # no idea what this is
layer1.grid ./ layer2.grid # elementwise division
./(layer1.grid, layer2.grid)

# Extend broadcast
# broadcast(f::Tf, As...) where Tf in Base.Broadcast
import Base.Broadcast: broadcast
function broadcast(f, layer::SimpleSDMLayer, As...)
    newlayer = copy(layer)
    newlayer.grid = broadcast(f, layer.grid, As...)
    return newlayer
end  

# Test extended broadcast
broadcast(+, layer1, 1)
layer1 + 1
layer1.grid .+ 1
@which broadcast(+, layer1.grid, 1)
broadcast(+, layer1, 1).grid
layer1 .+ 1
broadcast(+, layer1, layer2)
layer1 + layer2

##



