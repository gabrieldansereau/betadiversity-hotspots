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
    threshold = first(quantile(filter(!isnan, prediction[occ]), [0.05]))
    for i in eachindex(prediction.grid)
        prediction.grid[i] < threshold && (prediction.grid[i] = NaN)
    end
    return prediction
end

#=
function part2(prediction,occ)
    threshold = first(quantile(prediction[occ], [0.05]))
    for i in eachindex(prediction.grid)
        prediction.grid[i] < threshold && (prediction.grid[i] = NaN)
    end
    return prediction
end
warblers_test = warblers_occ[end-10:end]
@time predictions = pmap(x -> species_bclim(x, wc_vars_pred, train_vars = wc_vars_train), warblers_test)
problems = []
j = 0
for i in j+1:length(predictions)
    global j += 1
    part2(copy(predictions[i]), warblers_occ[end-10+i-1])
end
push!(problems, j)
# 2,3,4,5,7,8,9

j1 = 1
test1 = map(i -> part2(copy(predictions[i]), warblers_test[i]), j1:j1)
check1 = predictions[j1][warblers_test[j1]]
filtr1 = filter(!isnan, check1)
q1a = quantile(check1, [0.05])
q1b = quantile(filtr1, [0.05])

j2 = 2
test2 = map(i -> part2(copy(predictions[i]), warblers_test[i]), j2:j2)
check2 = predictions[j2][warblers_test[j2]]
filtr2 = filter(!isnan, check2)
q2a = quantile(check2, [0.05])
q2b = quantile(filtr2, [0.05])
=#
