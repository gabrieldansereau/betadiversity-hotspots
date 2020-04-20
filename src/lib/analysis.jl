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
    LCBDsets_raw = copy(LCBDsets)
    LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]
    LCBDsets = [LCBDsets..., LCBDsets_raw...]
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
