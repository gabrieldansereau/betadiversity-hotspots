import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# outcome = "sdm"

## Load presence-absence data for all species
@load "data/jld2/$(outcome)-distributions.jld2" distributions

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
using RCall
@rput Yobs
begin
    R"""
        library(vegan)
        Ytransf <- decostand(Yobs, "hel")
    """
end
@rget Ytransf

## Export results
# save_data = true
if (@defined save_data) && save_data == true
    @save "data/jld2/$(outcome)-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs
    @info "Data exported to file"
else
    @load "data/jld2/$(outcome)-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs
    @info "Data imported from file"
end

## Visualize patterns in Y
# Heatmap of Y
heat_y = heatmap(Yobs, title = "$(titlecase(outcome)) matrix Y (unsorted)",
                   ylabel = "Site number", xlabel = "Species number")
# Sort Y by rows & columns sums
rowsum = sum.(eachrow(Yobs))
colsum = sum.(eachcol(Yobs))
sortedrows = sortperm(rowsum)
sortedcols = sortperm(colsum, rev=true)
heat_sortrowcol = heatmap(Yobs[sortedrows, sortedcols], title = "$(titlecase(outcome)) matrix Y (sorted by row and column sums)",
                   ylabel = "Site number", xlabel = "Species number")

## Export results
# save_figures = true
if (@isdefined save_figures) && save_figure == true
    savefig(heat_sortrowcol, "fig/$(outcome)/02_$(outcome)_Y-rowcolsorted.png")
    @info "Figures saved"
else
    @info "Figures not saved"
end

#= # Funny looking smudge 😛
heatmap(sort(Yobs, dims=1, by=sum))
=#
