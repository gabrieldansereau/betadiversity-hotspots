# Convert vector values to quantile scores
function quantiles(vec::Array{Float64,1}; ignorenan::Bool=true)
    if ignorenan == true
        qfinder = ecdf(filter(!isnan, vec))
        qvec = replace(x -> !isnan(x) ? qfinder(x) : x, vec)
    else
        qfinder = ecdf(vec)
        qvec = qfinder(vec)
    end
    return qvec
end

# Convert & mutate vector values to quantile scores
function quantiles!(vec::Array{Float64,1}; ignorenan::Bool=true)
    if ignorenan == true
        qfinder = ecdf(filter(!isnan, vec))
        replace!(x -> !isnan(x) ? qfinder(x) : x, vec)
    else
        qfinder = ecdf(vec)
        replace!(x -> qfinder(x), vec)
    end
    return vec
end

# Convert SimpleSDMLayer values to quantile scores
function quantiles(layer::SimpleSDMLayer)
    qfinder = ecdf(filter(!isnan, layer.grid))
    qgrid = replace(x -> !isnan(x) ? qfinder(x) : x, layer.grid)
    qlayer = SimpleSDMResponse(qgrid, layer.left, layer.right, layer.bottom, layer.top)
    return qlayer
end

# Convert & mutate SimpleSDMLayer values to quantile scores
function quantiles!(layer::SimpleSDMLayer)
    qfinder = ecdf(filter(!isnan, layer.grid))
    replace!(x -> !isnan(x) ? qfinder(x) : x, layer.grid)
    return layer
end
