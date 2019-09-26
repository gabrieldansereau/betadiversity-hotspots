import SimpleSDMLayers: bioclim

function bioclim(occ::Union{GBIFRecords,DataFrame}, pred_vars::SimpleSDMLayer; train_vars::SimpleSDMLayer=pred_vars)
    observed_values = train_vars[occ]
    qfinder = ecdf(observed_values)
    lq = zeros(Float64, size(pred_vars))
    for i in eachindex(pred_vars.grid)
        if isnan(pred_vars.grid[i])
            local_quantile = NaN
        else
            local_quantile = qfinder(pred_vars.grid[i])
            if local_quantile > 0.5
                local_quantile = 1.0-local_quantile
            end
            local_quantile = 2.0 * local_quantile
        end
        lq[i] = local_quantile
    end
    prediction = SimpleSDMResponse(lq, pred_vars.left, pred_vars.right, pred_vars.bottom, pred_vars.top)
end
function species_bclim(occ, pred_vars; train_vars=pred_vars)
    predictions = [bioclim(occ, pred_vars[i], train_vars = train_vars[i]) for i in 1:length(pred_vars)];
    prediction = reduce(minimum, predictions);
    threshold = first(quantile(prediction[occ], [0.05]))
    for i in eachindex(prediction.grid)
        prediction.grid[i] < threshold && (prediction.grid[i] = NaN)
    end
    return prediction
end
