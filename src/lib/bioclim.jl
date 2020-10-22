## BIOCLIM model functions

# import SimpleSDMLayers: bioclim

# Small internal function for BIOCLIM. Ensures values are equal for both tails and between 0-1
_bcscore(x) = isnothing(x) ? nothing : (x > 0.5 ? 2*(1-x) : 2x)

# 1st part of BIOCLIM model for single environmental layer
function _bioclim_layer(occ::Union{G,D}, layer::T; training_layer::T=layer) where {T <: SimpleSDMLayer, D <: DataFrame, G <: GBIFRecords}
    # occ: occurences of a single species as a DataFrame with latitude and longitude columns
    # layer: layer of environmental variables used for prediction
    # training_layer: optional, training environmental layer, can use a different resolution

    # Get observed environmental values (training values)
    obs = training_layer[occ]
    filter!(!isnothing, obs)
    # Create empty layer for predictions
    pred = similar(layer)
    # Make predictions
    if length(unique(obs)) <= 1 # special case if 1 unique value only
        # If 1 unique value only, set 1.0 as prediction for all sites with this value
        pred.grid = replace(x -> x == obs[1] ? 1.0 : 0.0, layer.grid)
    else
        # Create ECDF function to extract quantile values
        qfinder = ecdf(obs)
        # Get prediction score
        qfgrid = replace(x -> isnothing(x) ? nothing : qfinder(x), layer.grid)
        pred.grid = _bcscore.(qfgrid)
    end
    # Restore nothings from original layer
    for idx in findall(isnothing, layer.grid)
        pred.grid[idx] = nothing
    end
    return pred
end

# Complete BIOCLIM model on all environmental layers, including 1st part function
function bioclim(occ, layers; training_layers=layers, binary=true, threshold::Float64=0.0)
    # occ: occurences of a single species as a DataFrame with latitude and longitude columns
    # layers: layers of environmental variables used for prediction
    # training_layers: optional, training environmental layers, can use a different resolution
    # binary: optional, convert result to binary presence-absence values
    # with_threshold: optional, apply threshold to remove lower predictions
    # threshold: optional, percentage of lower predictions to remove (e.g. 0.05)

    @assert 0.0 ≤ threshold < 1.0 "Threshold must be between 0.0 and 1.0"

    # Apply 1st part of BIOCLIM on each environmental variable
    predictions = [_bioclim_layer(occ, layers[i], training_layer = training_layers[i]) for i in eachindex(layers)];
    # Reduce to single layer with minimum values
    prediction = reduce(min, predictions);
    # Apply threshold (if specified)
    if threshold > 0.0
        # Select non-nothing values only
        no_nothing = filter(!isnothing, prediction[occ]);
        # Apply threshold only if there are non-nothing values
        if length(no_nothing) != 0
            # Get prediction value equivalent to threshold
            threshold_value = first(quantile(no_nothing, [threshold]));
            # Replace values smaller than threshold by nothing
            replace!(x -> !isnothing(x) || x < threshold_value ? nothing : x, prediction.grid);
        end
    end
    # Replace zeros by nothing
    replace!(x -> x == 0.0 ? nothing : x, prediction.grid);
    # Convert to binary presence-absence values (if requested, default is true)
    if binary == true
        replace!(x -> x > 0.0 ? 1.0 : x, prediction.grid)
    end
    return prediction
end
