## BIOCLIM model functions

# import SimpleSDMLayers: bioclim

# 1st part of BIOCLIM model for single environmental variables
function bioclim_singlevar(occ::Union{GBIFRecords,DataFrame}, pred_vars::SimpleSDMLayer; train_vars::SimpleSDMLayer=pred_vars)
    # occ: occurences of a single species as a DataFrame with latitude and longitude columns
    # pred_vars: environmental variables used for prediction
    # train_vars: optional, training environmental variables, can use a different resolution

    # Get observed environmental values (training values)
    observed_values = train_vars[occ]
    # Create ECDF function to extract quantile value
    qfinder = ecdf(observed_values)
    # Create empty array for local quantile values
    lq = zeros(Float64, size(pred_vars))
    # Loop for all sites (prediction values)
    for i in eachindex(pred_vars.grid)
        if isnan(pred_vars.grid[i])
            # Set value to NaN if env value is already NaN
            local_quantile = NaN
        else
            # Get quantile rank if not NaN
            local_quantile = qfinder(pred_vars.grid[i])
            # Replace values greater than 0.5 to set both tails equal
            if local_quantile > 0.5
                local_quantile = 1.0-local_quantile
            end
            # Scale back between 0 and 1
            local_quantile = 2.0 * local_quantile
        end
        # Collect quantile values
        lq[i] = local_quantile
    end
    # Convert to SimpleSDMLayer, with same coordinates as prediction layer
    prediction = SimpleSDMResponse(lq, pred_vars.left, pred_vars.right, pred_vars.bottom, pred_vars.top)
end

# Complete BIOCLIM model on all environmental variables, including 1st part function
function bioclim(occ, pred_vars; train_vars=pred_vars, binary=true, with_threshold=false, threshold::Float64=0.05)
    # occ: occurences of a single species as a DataFrame with latitude and longitude columns
    # pred_vars: environmental variables used for prediction
    # train_vars: optional, training environmental variables, can use a different resolution
    # binary: optional, convert result to binary presence-absence values
    # with_threshold: optional, apply threshold to remove lower predictions
    # threshold: optional, set threshold to use (default = 0.05)

    # Apply 1st part of BIOCLIM on each environmental variable
    predictions = [bioclim_singlevar(occ, pred_vars[i], train_vars = train_vars[i]) for i in 1:length(pred_vars)];
    # Reduce to single layer with minimum values
    prediction = reduce(minimum, predictions);
    # Apply threshold (if requested, default is false)
    if with_threshold == true
        # Select non-NaN values only
        no_nan = filter(!isnan, prediction[occ]);
        # Apply threshold only if there are non-NaN values
        if length(no_nan) != 0
            # Get prediction value equivalent to threshold
            threshold_value = first(quantile(no_nan, [threshold]));
            # Replace values smaller than threshold by NaN
            replace!(x -> x < threshold_value ? NaN : x, prediction.grid);
        end
    end
    # Replace zeros by NaN
    replace!(x -> x == 0.0 ? NaN : x, prediction.grid);
    # Convert to binary presence-absence values (if requested, default is true)
    if binary == true
        replace!(x -> x > 0.0 ? 1.0 : x, prediction.grid)
    end
    return prediction
end
