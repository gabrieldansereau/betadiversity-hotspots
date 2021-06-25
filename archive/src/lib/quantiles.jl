# Convert vector values to quantile scores
function quantiles(vec::Array; ignorenothing::Bool=true)
    if ignorenothing == true
        qfinder = ecdf(filter(!isnothing, vec))
        qvec = replace(x -> !isnothing(x) ? qfinder(x) : x, vec)
    else
        qfinder = ecdf(vec)
        qvec = qfinder(vec)
    end
    return qvec
end

# Convert & mutate vector values to quantile scores
function quantiles!(vec::Array; ignorenothing::Bool=true)
    if ignorenothing == true
        qfinder = ecdf(filter(!isnothing, vec))
        replace!(x -> !isnothing(x) ? qfinder(x) : x, vec)
    else
        qfinder = ecdf(vec)
        replace!(x -> qfinder(x), vec)
    end
    return vec
end

# Convert SimpleSDMLayer values to quantile scores
function quantiles(layer::SimpleSDMLayer)
    qfinder = ecdf(Float32.(filter(!isnothing, layer.grid)))
    qgrid = replace(x -> isnothing(x) ? x : qfinder(x), layer.grid)
    qlayer = SimpleSDMResponse(qgrid, layer)
    return qlayer
end

# Convert & mutate SimpleSDMLayer values to quantile scores
function quantiles!(layer::SimpleSDMLayer)
    qfinder = ecdf(Float32.(filter(!isnothing, layer.grid)))
    replace!(x -> isnothing(x) ? x : qfinder(x), layer.grid)
    return layer
end
