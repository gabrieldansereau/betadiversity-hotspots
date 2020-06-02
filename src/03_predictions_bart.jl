import Pkg
Pkg.activate(".")
using RCall
begin
    R"""
    source(file.path("src", "required.R"))
    """
end
using Distributed
@time include("required.jl")

## Conditional arguments
# save_data = true

## Perform BARTs
#=
@time begin
    R"""
    create_models <- FALSE
    save_models <- FALSE

    source(here("src", "02_training_bart.R"))
    """
end
@rget pred_df lower_df upper_df pres_df results
=#

# Load predictions from CSV
results  = CSV.read(joinpath("data", "proc", "bart_summaries.csv"))
varimps  = CSV.read(joinpath("data", "proc", "bart_varimps.csv"))
pred_df  = CSV.read(joinpath("data", "proc", "bart_predictions_prob.csv"), missingstrings = ["NA"])
lower_df = CSV.read(joinpath("data", "proc", "bart_predictions_lower.csv"), missingstrings = ["NA"])
upper_df = CSV.read(joinpath("data", "proc", "bart_predictions_upper.csv"), missingstrings = ["NA"])
pres_df  = CSV.read(joinpath("data", "proc", "bart_predictions_pres.csv"), missingstrings = ["NA"])


## Create Y matrices

# Get matrix Y
Y = replace(Array(pres_df), missing => NaN)
Yprob  = replace(Array(pred_df),  missing => NaN)
Ylower = replace(Array(lower_df), missing => NaN)
Yupper = replace(Array(upper_df), missing => NaN)
# Set values to NaN if no species present
inds_zeros = findall(map(x -> all(iszero.(x)), eachrow(Y)))
Y[inds_zeros,:] .= NaN

## Create distributions

# Load raw distributions (for grid size)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
raw_distributions = distributions
# Cut to Quebec coordinates (optional)
# coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
# raw_distributions = [d[coords_qc] for d in raw_distributions]
# Get layer dimensions & limits
dims = size(raw_distributions[1].grid)
lims = (left = raw_distributions[1].left, right = raw_distributions[1].right,
        bottom = raw_distributions[1].bottom, top = raw_distributions[1].top)

# Create distribution layers
layers = []
for Y in (Y, Yprob, Ylower, Yupper)
    Ydistrib = replace(Y, 0.0 => NaN)
    Ygrids = [Ydistrib[:, col] for col in 1:size(Ydistrib,2)]
    Ygrids = reshape.(Ygrids, dims...)
    distributions = SimpleSDMResponse.(Ygrids, lims...)
    push!(layers, distributions)
end
distributions, prob_distrib, lower_distrib, upper_distrib = layers;
distributions

## Export results
# save_data = true
if (@isdefined save_data) && save_data == true
    @save joinpath("data", "jld2", "bart-distributions.jld2") distributions prob_distrib lower_distrib upper_distrib
    _zip_jld2(joinpath("data", "jld2", "bart-distributions.zip"),
              joinpath("data", "jld2", "bart-distributions.jld2"))
    touch(joinpath("data", "jld2", "bart-distributions.jld2"))
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])
