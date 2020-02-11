import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# Load old raw
@load "/home/gdansereau/pcloud/raw-distributions.jld2" distributions spenames speindex
@load "/home/gdansereau/pcloud/raw-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
oldraw = (distributions = distributions,
          spenames = spenames,
          speindex = speindex,
          Y = Y,
          Yobs = Yobs,
          Ytransf = Ytransf,
          inds_obs = inds_obs,
          inds_notobs = inds_notobs,
          pres_counts = pres_counts)

# Load old SDM
@load "/home/gdansereau/pcloud/sdm-distributions.jld2" distributions spenames speindex
@load "/home/gdansereau/pcloud/sdm-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
oldsdm = (distributions = distributions,
          spenames = spenames,
          speindex = speindex,
          Y = Y,
          Yobs = Yobs,
          Ytransf = Ytransf,
          inds_obs = inds_obs,
          inds_notobs = inds_notobs,
          pres_counts = pres_counts)

# Load new raw
@load "/home/gdansereau/github/BioClim/data/jld2/raw-distributions.jld2" distributions spenames speindex
@load "/home/gdansereau/github/BioClim/data/jld2/raw-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
newraw = (distributions = distributions,
          spenames = spenames,
          speindex = speindex,
          Y = Y,
          Yobs = Yobs,
          Ytransf = Ytransf,
          inds_obs = inds_obs,
          inds_notobs = inds_notobs,
          pres_counts = pres_counts)

# Load new SDM
@load "/home/gdansereau/github/BioClim/data/jld2/sdm-distributions.jld2" distributions spenames speindex
@load "/home/gdansereau/github/BioClim/data/jld2/sdm-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
newsdm = (distributions = distributions,
          spenames = spenames,
          speindex = speindex,
          Y = Y,
          Yobs = Yobs,
          Ytransf = Ytransf,
          inds_obs = inds_obs,
          inds_notobs = inds_notobs,
          pres_counts = pres_counts)

distributions = nothing
spenames = nothing
speindex = nothing
Y = nothing
Yobs = nothing
Ytransf = nothing
inds_obs = nothing
inds_notobs = nothing
pres_counts = nothing

##

countsdf = DataFrame(oldraw = oldraw.pres_counts,
                     newraw = newraw.pres_counts,
                     diffraw = newraw.pres_counts .- oldraw.pres_counts,
                     oldsdm = oldsdm.pres_counts,
                     newsdm = newsdm.pres_counts,
                     diffsdm = newsdm.pres_counts .- oldsdm.pres_counts,)

show(countsdf, allrows=true)
