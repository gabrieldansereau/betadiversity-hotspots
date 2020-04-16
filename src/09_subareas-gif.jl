import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

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
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

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
include(joinpath("lib", "beta-div.jl"))
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

## Combine figures
function plot_lcbd_richness(richness, lcbd; title = "", kw...)
  p1 = plot(richness, c = :viridis, title = "Richness", colorbar_title = "Number of species")
  p2 = plot(lcbd, c = :viridis, title = "LCBD", colorbar_title = "Relative LCBD score")
  p3 = plot(quantiles(lcbd), c = :viridis, title = "LCBD quantiles", colorbar_title = "Quantile rank")
  p4 = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
            xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites")
  if title != ""
    l = @layout [t{.01h}; grid(2,2)]
    ptitle = plot(annotation = (0.5, 0.5, "$title"), framestyle = :none)
    p = plot(ptitle, p1, p2, p4, p3, layout = l; kw...)
  else
    p = plot(p1, p2, p4, p3; kw...)
  end
  return p
end

#### Repeat for different subareas
function plot_subareas(coords, initial_distributions; transform = true, kw...)
  distributions = [d[coords] for d in initial_distributions]
  Y = calculate_Ymatrix(distributions)
  richness = calculate_richness(Y.Y, Y.inds_notobs, distributions)
  lcbd = calculate_lcbd(Y.Yobs, Y.Ytransf, Y.inds_obs, distributions)
  if transform
    p = plot_lcbd_richness(richness, lcbd[2]; kw...)
  else
    p = plot_lcbd_richness(richness, lcbd[1]; kw...)
  end
end

left = -71.0; right = -64.0; bottom = 47.5; top = 50.0
coords_subarea = (left = left, right = right, bottom = bottom, top = top)
p = plot_subareas(coords_subarea, distributions)
p = plot_subareas(coords_subarea, distributions; formatter = f -> "$(round(f, digits = 1))")

subarea_plots = []

@time while left > -145.0 && bottom > 20.0
  global left -= 1.0
  global bottom -= 0.5;
  coords_subarea = (left = left, right = right, bottom = bottom, top = top)
  p = plot_subareas(coords_subarea, distributions, formatter = f -> "$(round(f, digits = 1))")
  push!(subarea_plots, p)
end
subarea_plots

anim = @animate for p in subarea_plots
    plot(p)
end
gif(anim, joinpath("fig", outcome, "09_subareas.gif"), fps = 3)
