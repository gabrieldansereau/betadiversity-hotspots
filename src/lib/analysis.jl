## Y matrix
function Ymatrix(distributions; transform=false, observed=false)
    ## Create matrix Y (site-by-species community data table)
    # Get distributions as vectors
    distributions_vec = [vec(d.grid) for d in distributions]
    # Create matrix Y by combining distribution vectors
    Y = hcat(distributions_vec...)

    # Get indices of sites with observations
    inds_obs = _indsobs(Y)
    # Create matrix Yobs with observed sites only, replace nothings by zeros
    Yobs = _Yobs(Y, inds_obs)
    # Apply Hellinger transformation (using vegan in R)
    if transform
        Yobs = _Ytransf(Yobs)
    end
    # Replace values in original matrix
    Y[inds_obs, :] = Yobs

    if observed
        return Yobs
    else
        return Y
    end
end

function _indsobs(Y)
    # Verify if sites have observations
    sites_obs = [any(y -> !isnothing(y) && y > 0, yrow) for yrow in eachrow(Y)]
    # Get indices of sites with observations
    inds_obs = findall(sites_obs)
    return inds_obs
end

function _indsnotobs(Y)
    # Verify if sites have observations
    sites_obs = [any(y -> !isnothing(y) && y > 0, yrow) for yrow in eachrow(Y)]
    # Get indices of sites without observations
    inds_notobs = findall(.!sites_obs)
    return inds_notobs
end

function _Yobs(Y, inds_obs)
    # Create matrix Yobs with observed sites only
    Yobs = Y[inds_obs, :]
    # Replace nothings by zeros for observed sites (~true absences)
    replace!(Yobs, nothing => 0.0)
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
function calculate_richness(Y, layer)
    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)
    # Create empty empty vector
    sums = fill(nothing, size(Y, 1)) |> Array{Union{Nothing,Float32}}
    # Get number of species per observed site
    sums_obs = map(sum, eachrow(Yobs))
    sums[inds_obs] = sums_obs
    # Reshape to grid format
    sums = reshape(sums, size(layer)) |> Array
    ## Create SimpleSDMLayer
    return richness = SimpleSDMResponse(sums, layer)
end

function calculate_gamma(Y)
    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)
    sp_counts = sum(Yobs; dims=1)
    gamma = sum(sp_counts .> 0)
    return gamma
end

## LCBD
# Load functions
function calculate_lcbd(Y, layer; transform=true, relative=true)
    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)

    # Apply hellinger transformation
    if transform
        Yobs = _Ytransf(Yobs)
    end
    # Compute beta diversity statistics
    BDstats = betadiv(Yobs)

    # Extract LCBD values
    LCBDvals = BDstats.LCBDi
    # Scale LCBDi values to maximum value
    if relative
        LCBDvals = LCBDvals ./ maximum(LCBDvals)
    end

    # Create empty grid
    LCBDgrid = fill(nothing, size(layer)) |> Array{Union{Nothing,Float32}}
    # Fill-in grid
    LCBDgrid[inds_obs] = LCBDvals
    # Create SimpleSDMLayer with LCBD values
    LCBDlayer = SimpleSDMResponse(LCBDgrid, layer)
    return LCBDlayer
end

function calculate_BDtotal(Y; transform=true)
    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)

    # Apply hellinger transformation
    if transform
        Yobs = _Ytransf(Yobs)
    end
    # Compute beta diversity statistics
    BDstats = betadiv(Yobs)

    # Extract LCBD values
    BDtotal = BDstats.BDtotal

    return BDtotal
end
