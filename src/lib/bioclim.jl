import SimpleSDMLayers: bioclim

function bioclim(p::SimpleSDMLayer, r::Union{GBIFRecords,DataFrame})
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
    prediction = SimpleSDMResponse(lq, p.left, p.right, p.bottom, p.top)
end
function bioclim_training(t::SimpleSDMLayer, r::Union{GBIFRecords,DataFrame})
    observed_values = t[r]
    qfinder = ecdf(observed_values)
    return qfinder
end
function bioclim_prediction(p::SimpleSDMLayer, r::Union{GBIFRecords,DataFrame}, qfinder::ECDF{Array{Float64,1}})
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
    prediction = SimpleSDMResponse(lq, p.left, p.right, p.bottom, p.top)
end
function species_bclim(occ, vars)
    predictions = [bioclim(vars[i], occ) for i in 1:length(vars)];
    prediction = reduce(minimum, predictions);
    for i in eachindex(prediction.grid)
        if prediction.grid[i] == 0.0
            prediction.grid[i] = NaN
        end
    end
    return prediction
end
