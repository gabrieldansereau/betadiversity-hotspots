## Y matrix function

"""
    Ymatrix(distributions::Vector{<:SimpleSDMLayer}; transform::Bool=false, observed::Bool=false)

Creates a matrix Y (a site-by-species community data table) from a vector of
layers with species distributions. Sites (rows) are ordered as in the layer
grids and species (columns) are ordered as in the vector of layers.

In the matrix Y, values for sites where at least one species is present are
either 0.0 for absence and 1.0 for presence. On the other hand, values for sites
where no species are present are all `nothing`. These sites can be removed with
`observed=true` (default is `false`).

Setting `transform=true` (default is `false`) will apply the Hellinger
transformation to the Y matrix using the function `decostand` from the package
`vegan` in *R*. However, as this step is only needed to compute LCBD scores, it
is recommended to keep the default and use the `transform` argument from the
`lcbd` function instead.
"""
function Ymatrix(
    distributions::Vector{<:SimpleSDMLayer}; transform::Bool=false, observed::Bool=false
)
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

## Richness functions

"""
    richness(Y::Matrix, layer::SimpleSDMLayer)

Computes the species richness for the sites in the community matrix `Y` and
returns a `SimpleSDMResponse` with the same dimensions and coordinates and
dimensions as `layer`.
"""
function richness(Y::Matrix, layer::SimpleSDMLayer)
    @assert size(Y, 1) == length(layer.grid) "Y and layer must have the same number of sites"
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
    return SimpleSDMResponse(sums, layer)
end

"""
    gamma(Y::Matrix)

Returns the gamma diversity in the region covered by the matrix `Y`, i.e. the
number of species with a least one presence in `Y`. Note that this will differ
from the number of columns when some species have no observations in the region.
"""
function gamma(Y::Matrix)
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)
    sp_counts = sum(Yobs; dims=1)
    return sum(sp_counts .> 0)
end

## Beta diversity function

"""
    lcbd(Y::Matrix, layer::SimpleSDMLayer; transform::Bool=true, relative::Bool=true)

Computes the LCBD (local contributions to beta diversity) scores for the sites
in the matrix `Y` and returns a `SimpleSDMResponse` with the same dimensions and
coordinates as `layer`.

Setting `transform=true` (the default) will apply the Hellinger transformation
to the Y matrix using the function `decostand` from the package `vegan` in *R*,
which is appropriate for presence-absence data.

Setting `relative=true` (default is `false`) will rescale the LCBD scores as
relative to the maximum score (whose value will be one).
"""
function lcbd(Y::Matrix, layer::SimpleSDMLayer; transform::Bool=true, relative::Bool=false)
    @assert size(Y, 1) == length(layer.grid) "Y and layer must have the same number of sites"

    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)

    # Apply hellinger transformation (if requested)
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

"""
    beta_total(Y::Matrix; transform::Bool=true)

Returns the total beta diversity from the matrix `Y`, i.e. the unbiased &
comparable estimator of variance in `Y`.
"""
function beta_total(Y::Matrix; transform::Bool=true)
    # Create necessary Y elements
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)

    # Apply hellinger transformation
    if transform
        Yobs = _Ytransf(Yobs)
    end
    # Compute beta diversity statistics
    BDstats = betadiv(Yobs)

    return BDstats.BDtotal
end
