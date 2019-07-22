function bioclim(p::SDMLayer, r::Union{GBIFRecords,DataFrame})
    observed_values = p[r]
    qfinder = ecdf(observed_values)
    lq = zeros(Float64, size(p))
    for i in eachindex(p.grid)
        if isnan(p.grid[i])
            local_quantile = NaN
        else
            local_quantile = qfinder(p.grid[i])
            if local_quantile > 0.5
                local_quantile = 1.0-local_quantile
            end
            local_quantile = 2.0 * local_quantile
        end
        lq[i] = local_quantile
    end
    prediction = SDMLayer(lq, p.left, p.right, p.bottom, p.top)
end
function bioclim_training(t::SDMLayer, r::Union{GBIFRecords,DataFrame})
    observed_values = t[r]
    qfinder = ecdf(observed_values)
    return qfinder
end
function bioclim_prediction(p::SDMLayer, r::Union{GBIFRecords,DataFrame}, qfinder::ECDF{Array{Float64,1}}) 
    lq = zeros(Float64, size(p))
    for i in eachindex(p.grid)
        if isnan(p.grid[i])
            local_quantile = NaN
        else
            local_quantile = qfinder(p.grid[i])
            if local_quantile > 0.5
                local_quantile = 1.0-local_quantile
            end
            local_quantile = 2.0 * local_quantile
        end
        lq[i] = local_quantile
    end
    prediction = SDMLayer(lq, p.left, p.right, p.bottom, p.top)
end
