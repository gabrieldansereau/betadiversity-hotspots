using Distributed
using JLD2
@time include("../required.jl")

@time begin
    # Load data from CSV files
    df = CSV.read("data/proc/ebd_warblers_prep.csv", header=true, delim="\t")
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]
end


@time pres_abs = @showprogress pmap(x -> presence_absence(x, wc_vars[1], binary = false), warblers_occ)

pres_abs

## Create matrix Y (site-by-species community data table)
begin
    # Get dimensions
    nsites = prod(size(pres_abs[1]))
    nspecies = length(pres_abs)
    # Create Y
    Y = zeros(Float64, (nsites, nspecies))
    # Fill Y with community predictions
    @progress for gc in 1:nsites # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], pres_abs)
        # Fill Y with binary values
        global Y[gc,:] = R
    end
end

## Create matrix Ypred (only sites with observations)
# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)
inds_notpred = findall(.!sites_pred)
# Select sites with predictions only
Ypred = Y[inds_pred,:]

mean_count_per_obs = map(x -> mean(filter(!iszero, x)), eachcol(Ypred))
nb_sites_with_obs = map(x -> length(filter(!iszero, x)), eachcol(Ypred))
species = unique(df.species)
burstiness = mean_count_per_obs./nb_sites_with_obs

test = DataFrame(species = species, mean_count_per_obs = mean_count_per_obs, nb_sites_with_obs = nb_sites_with_obs)
test.burstiness = test.mean_count_per_obs./test.nb_sites_with_obs
test
show(sort(test, :burstiness, rev=true), allrows=true)
show(sort(test, :mean_count_per_obs, rev=true), allrows=true)
