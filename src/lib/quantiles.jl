# Function to extract quantile scores of a vector
function quantiles(vec::Array{Float64,1})
    ecdf_finder = ecdf(vec)
    return ecdf_finder(vec)
end

# Function to extract quantile scores of a SimpleSDMLayer
function quantiles(layer::SimpleSDMLayer)
    ecdf_finder = ecdf(filter(!isnan, layer.grid))
    qgrid = replace(x -> !isnan(x) ? ecdf_finder(x) : x, layer.grid)
    qlayer = SimpleSDMResponse(qgrid, layer.left, layer.right, layer.bottom, layer.top)
    return qlayer
end
