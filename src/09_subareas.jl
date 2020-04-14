import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# Make sure "outcome" is defined
outcome = "rf"
if !(@isdefined outcome)
  @warn "'outcome' not defined, must be either 'raw', 'sdm' or 'rf'"
elseif !(outcome in ["raw", "sdm", "rf"])
  @warn "'outcome' invalid, must be either 'raw', 'sdm' or 'rf'"
else
  @info "'outcome' currently set to '$(outcome)'"
end

## Load distribution data for all species
@load "data/jld2/$(outcome)-distributions.jld2" rf_distributions
distributions = rf_distributions

## Extract subareas
# Northeast subarea
coords_NE = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
distributions_NE = [d[coords_NE] for d in distributions]
# Southwest subarea
coords_SW = (left = -120.0, right = -100.0, bottom = 30.0, top = 40.0)
distributions_SW = [d[coords_SW] for d in distributions]

## Get Ymatrices
using RCall
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
NE = calculate_Ymatrix(distributions_NE)
SW = calculate_Ymatrix(distributions_SW)


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
richness_NE = calculate_richness(NE.Y, NE.inds_notobs, distributions_NE)
richness_SW = calculate_richness(SW.Y, SW.inds_notobs, distributions_SW)

# Visualize
plot(richness_NE, c = :viridis)
plot(richness_SW, c = :viridis)

## LCBD
# Load functions
include("lib/beta-div.jl")
function calculate_lcbd(Yobs, Ytransf, inds_obs, distributions)
  ## Compute beta diversity statistics
  # Compute BD statistics on distribution data
  resBDobs = BD(Yobs)
  # Compute BD statistics on transformed data
  resBDtransf = BD(Ytransf)

  # Extract LCBD values
  resBD = [resBDobs, resBDtransf]
  LCBDsets = [res.LCBDi for res in resBD]
  # Scale LCBDi values to maximum value
  LCBDsets_raw = copy(LCBDsets)
  LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]
  LCBDsets = [LCBDsets..., LCBDsets_raw...]

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
lcbd_NE = calculate_lcbd(NE.Yobs, NE.Ytransf, NE.inds_obs, distributions_NE)
lcbd_SW = calculate_lcbd(SW.Yobs, SW.Ytransf, SW.inds_obs, distributions_SW)

# Visualize
plot(lcbd_NE[1], c = :viridis)
plot(lcbd_NE[2], c = :viridis)
plot(lcbd_SW[1], c = :viridis)
plot(lcbd_SW[2], c = :viridis)
# Quantiles
plot(quantiles(lcbd_NE[1]), c = :viridis)
plot(quantiles(lcbd_NE[2]), c = :viridis)
plot(quantiles(lcbd_SW[1]), c = :viridis)
plot(quantiles(lcbd_SW[2]), c = :viridis)

## Relationship
rel_NE = histogram2d(richness_NE, lcbd_NE[1], c = :viridis, bins = 40,
          xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
rel_tr_NE = histogram2d(richness_NE, lcbd_NE[2], c = :viridis, bins = 40,
              xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
rel_SW = histogram2d(richness_SW, lcbd_SW[1], c = :viridis, bins = 40,
          xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
rel_tr_SW = histogram2d(richness_SW, lcbd_SW[2], c = :viridis, bins = 40,
              xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
