# import SimpleSDMLayers: bioclim

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
function species_bclim(occ, pred_vars; train_vars=pred_vars, binary=true, with_threshold=false)
    predictions = [bioclim(occ, pred_vars[i], train_vars = train_vars[i]) for i in 1:length(pred_vars)];
    prediction = reduce(minimum, predictions);
    if with_threshold == true
        no_nan = filter(!isnan, prediction[occ]);
        if length(no_nan) != 0
            threshold = first(quantile(no_nan, [0.05]));
            replace!(x -> x < threshold ? NaN : x, prediction.grid);
        end
    end
    replace!(x -> x == 0.0 ? NaN : x, prediction.grid);
    if binary == true
        replace!(x -> x > 0.0 ? 1.0 : x, prediction.grid)
    end
    return prediction
end
