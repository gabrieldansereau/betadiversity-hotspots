## Y matrix
function calculate_Y(distributions; transform = false)
    ## Create matrix Y (site-by-species community data table)
    # Get distributions as vectors
    distributions_vec = [vec(d.grid) for d in distributions];
    # Create matrix Y by combining distribution vectors
    Y = hcat(distributions_vec...);

    # Get indices of sites with observations
    inds_obs = _indsobs(Y)
    # Create matrix Yobs with observed sites only, replace nothings by zeros
    Yobs = _Yobs(Y, inds_obs)
    # Apply Hellinger transformation (using vegan in R)
    if transform
        Yobs = _Ytransf(Yobs)
    end
    # Replace values in original matrix
    Y[inds_obs,:] = Yobs;

    return Y
end

function _indsobs(Y)
    # Verify if sites have observations
    sites_obs = [any(y -> !isnothing(y) && y > 0, yrow) for yrow in eachrow(Y)];
    # Get indices of sites with observations
    inds_obs = findall(sites_obs);
    return inds_obs
end

function _indsnotobs(Y)
    # Verify if sites have observations
    sites_obs = [any(y -> !isnothing(y) && y > 0, yrow) for yrow in eachrow(Y)];
    # Get indices of sites without observations
    inds_notobs = findall(.!sites_obs);
    return inds_notobs
end

function _Yobs(Y, inds_obs)
    # Create matrix Yobs with observed sites only
    Yobs = Y[inds_obs,:];
    # Replace NaNs by zeros for observed sites (~true absences)
    replace!(Yobs, nothing => 0.0);
    return Yobs
end
_Yobs(Y) = _Yobs(Y, _indsobs(Y))

function _Ytransf(Yobs)
    # Remove type Nothing from Array (weird effects with RCall)
    Yobs = Array{Float32}(Yobs)
    # Apply Hellinger transformation (using vegan in R)
    @rput Yobs
    begin
        R"""
            Ytransf <- vegan::decostand(Yobs, "hel")
        """
    end
    @rget Ytransf
    return Ytransf
end

## Richness
function calculate_richness(Y, inds_notobs, layer)
    ## Get number of species per site
    sums = map(x -> Float64(sum(x)), eachrow(Y))
    # Add NaN for non predicted sites
    sums[inds_notobs] .= NaN
    # Reshape to grid format
    sums = reshape(sums, size(layer))
    ## Create SimpleSDMLayer
    richness = SimpleSDMResponse(sums, layer.left, layer.right, layer.bottom, layer.top)
end
calculate_richness(Y, layer) = calculate_richness(Y, _indsnotobs(Y), layer)

## LCBD
# Load functions
function calculate_lcbd(Yobs, inds_obs, layer; transform = true, relative = true)
    # Apply hellinger transformation
    if transform
        Yobs = _Ytransf(Yobs)
    end
    # Compute beta diversity statistics
    BDstats = BD(Yobs)

    # Extract LCBD values
    LCBDvals = BDstats.LCBDi
    # Scale LCBDi values to maximum value
    if relative
        LCBDvals = LCBDvals ./ maximum(LCBDvals)
    end

    # Create empty grid
    LCBDgrid = fill(NaN, size(layer))
    # Fill-in grid
    LCBDgrid[inds_obs] = LCBDvals
    # Create SimpleSDMLayer with LCBD values
    LCBDlayer = SimpleSDMResponse(LCBDgrid, layer.left, layer.right,
                                    layer.bottom, layer.top)
    return LCBDlayer
end

function calculate_lcbd(Y, layer; kw...)
    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)
    # Compute LCBD indices
    calculate_lcbd(Yobs, inds_obs, layer; kw...)
end
