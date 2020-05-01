## Get Ymatrices
function calculate_Ymatrix(distributions)
  ## Create matrix Y (site-by-species community data table)
  # Get distributions as vectors
  distributions_vec = [vec(d.grid) for d in distributions];
  # Create matrix Y by combining distribution vectors
  Y = hcat(distributions_vec...);

  # Verify if sites have observations
  sites_obs = [any(y .> 0.0) for y in eachrow(Y)];
  # Get indices of sites with observations
  inds_obs = findall(sites_obs);
  # Get indices of sites without observations
  inds_notobs = findall(.!sites_obs);

  # Create matrix Yobs with observed sites only
  Yobs = Y[inds_obs,:];
  # Replace NaNs by zeros for observed sites (~true absences)
  replace!(Yobs, NaN => 0.0);
  # Replace NaNs in original matrix Y too
  Y[inds_obs,:] = Yobs;

  ## Apply Hellinger transformation (using vegan in R)
  @rput Yobs
  begin
      R"""
          library(vegan)
          Ytransf <- decostand(Yobs, "hel")
      """
  end
  @rget Ytransf

  output = (Y = Y, Yobs = Yobs, Ytransf = Ytransf,
            inds_obs = inds_obs, inds_notobs = inds_notobs)
  return output
end

function calculate_Ymatrix(distributions; transform = false)
  ## Create matrix Y (site-by-species community data table)
  # Get distributions as vectors
  distributions_vec = [vec(d.grid) for d in distributions];
  # Create matrix Y by combining distribution vectors
  Y = hcat(distributions_vec...);

  # Get indices of sites with observations
  inds_obs = _indsobs(Y)
  # Create matrix Yobs with observed sites only, replace NaNs by zeros
  Yobs = _Yobs(Y)
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
  sites_obs = [any(y .> 0.0) for y in eachrow(Y)];
  # Get indices of sites with observations
  inds_obs = findall(sites_obs);
  return inds_obs
end

function _indsnotobs(Y)
  # Verify if sites have observations
  sites_obs = [any(y .> 0.0) for y in eachrow(Y)];
  # Get indices of sites without observations
  inds_notobs = findall(.!sites_obs);
  return inds_notobs
end

function _Yobs(Y)
  inds_obs = _indsobs(Y)
  # Create matrix Yobs with observed sites only
  Yobs = Y[inds_obs,:];
  # Replace NaNs by zeros for observed sites (~true absences)
  replace!(Yobs, NaN => 0.0);
  return Yobs
end

function _Ytransf(Yobs)
  ## Apply Hellinger transformation (using vegan in R)
  @rput Yobs
  begin
      R"""
          library(vegan)
          Ytransf <- decostand(Yobs, "hel")
      """
  end
  @rget Ytransf
  return Ytransf
end

## Richness
function calculate_richness(Y, inds_notobs, distributions)
  ## Get number of species per site
  sums = map(x -> Float64(sum(x)), eachrow(Y))
  # Add NaN for non predicted sites
  sums[inds_notobs] .= NaN
  # Reshape to grid format
  sums = reshape(sums, size(distributions[1]))
  ## Create SimpleSDMLayer
  richness = SimpleSDMResponse(sums, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)
end

calculate_richness(Y, distributions) = calculate_richness(Y, _indsnotobs(Y), distributions)

## LCBD
# Load functions
function calculate_lcbd(Yobs, Ytransf, inds_obs, distributions; relative = true)
  ## Compute beta diversity statistics
  # Compute BD statistics on distribution data
  resBDobs = BD(Yobs)
  # Compute BD statistics on transformed data
  resBDtransf = BD(Ytransf)

  # Extract LCBD values
  resBD = [resBDobs, resBDtransf]
  LCBDsets = [res.LCBDi for res in resBD]
  # Scale LCBDi values to maximum value
  if relative
    LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]
  end

  ## Arrange LCBD values as grid
  # Create empty grids
  LCBDgrids = [fill(NaN, size(distributions[1])) for LCBDi in LCBDsets]
  # Fill in grids
  [LCBDgrids[i][inds_obs] = LCBDsets[i] for i in 1:length(LCBDgrids)]
  # Create SimpleSDMLayer with LCBD values
  LCBD = SimpleSDMResponse.(LCBDgrids, distributions[1].left, distributions[1].right,
                            distributions[1].bottom, distributions[1].top)
  return LCBD
end

function calculate_lcbd(Y, distributions; kw...)
  # Create necessary Y elements
  Yobs = _Yobs(Y)
  Ytransf = _Ytransf(Yobs)
  inds_obs = _indsobs(Y)
  # Compute LCBD indices
  calculate_lcbd(Yobs, Ytransf, inds_obs, distributions; kw...)
end
