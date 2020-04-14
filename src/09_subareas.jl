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

coords_subarea1 = (left = -80.0, right = -60.0, bottom = 40.0, top = 50.0)
coords_subarea2 = (left = -120.0, right = -100.0, bottom = 30.0, top = 40.0)

distributions_sa1 = [d[coords_subarea1] for d in distributions]
distributions_sa2 = [d[coords_subarea2] for d in distributions]

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

sa1 = calculate_Ymatrix(distributions_sa1)
sa2 = calculate_Ymatrix(distributions_sa2)


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

richness_sa1 = calculate_richness(sa1.Y, sa1.inds_notobs, distributions_sa1)
richness_sa2 = calculate_richness(sa2.Y, sa2.inds_notobs, distributions_sa2)

plotSDM(richness_sa1, c = :viridis)
plotSDM(richness_sa2, c = :viridis)

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
lcbd_sa1 = calculate_lcbd(sa1.Yobs, sa1.Ytransf, sa1.inds_obs, distributions_sa1)
lcbd_sa2 = calculate_lcbd(sa2.Yobs, sa2.Ytransf, sa2.inds_obs, distributions_sa2)

plotSDM(lcbd_sa1[1], c = :viridis)
plotSDM(lcbd_sa1[2], c = :viridis)
plotSDM(lcbd_sa2[1], c = :viridis)
plotSDM(lcbd_sa2[2], c = :viridis)

plotSDM(quantiles(lcbd_sa1[1]), c = :viridis)
plotSDM(quantiles(lcbd_sa1[2]), c = :viridis)
plotSDM(quantiles(lcbd_sa2[1]), c = :viridis)
plotSDM(quantiles(lcbd_sa2[2]), c = :viridis)

histogram2d(richness_sa1, lcbd_sa1[1], c = :viridis, bins = 40,
            xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
histogram2d(richness_sa1, lcbd_sa1[2], c = :viridis, bins = 40,
            xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
histogram2d(richness_sa2, lcbd_sa2[1], c = :viridis, bins = 40, 
            xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
histogram2d(richness_sa2, lcbd_sa2[2], c = :viridis, bins = 40,
            xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
