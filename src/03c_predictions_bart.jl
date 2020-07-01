#### 03c - BART predictions ####
import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

## Conditional arguments
# save_data = true


## Load predictions from CSV
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
    # Distributions only
    jld_path = joinpath("data", "jld2", "bart-distributions.jld2")
    @save jld_path distributions
    _zip_jld2(replace(jld_path, ".jld2" => ".zip"), jld_path)
    touch(jld_path)
    # Extras
    jld_path = joinpath("data", "jld2", "bart-distributions_xtras.jld2")
    @save jld_path prob_distrib lower_distrib upper_distrib
    _zip_jld2(replace(jldpath, ".jld2" => ".zip"), jld_path)
    touch(jld_path)
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])

## Map uncertainty

# Get uncertainty per species
uncertainty = upper_distrib .- lower_distrib
plotSDM2(uncertainty[1], c = :viridis)

# Uncertainty sum
uncertainty_sum = sum(uncertainty)
plotSDM2(uncertainty_sum, c = :viridis)
# Uncertainty mean
uncertainty_mean = mean(uncertainty)
plotSDM2(uncertainty_mean, c = :viridis)
histogram(uncertainty_mean)

# Plot uncertainty & histogram
function uncertainty_plot(layer; title = "")
    unc_map = plotSDM2(layer,
                        c = :viridis, clim = extrema(layer),
                        title = "Uncertainty",
                        colorbar_title = "Uncertainty mean")
    unc_hist = histogram(layer,
                          legend = :none,
                          # ylim = maxrange, # xlabel = "Difference",
                          title = "Distribution of uncertainty values",
                          orientation = :horizontal)
    unc_title = plot(annotation = (0.5, 0.5, "$(title)"), framestyle = :none)
    l = @layout [t{0.01h}; a{0.6w} b{0.38w}]
    unc_plot = plot(unc_title, unc_map, unc_hist,
                     size = (800, 400), layout = l)
    return unc_plot
end
unc_plot = uncertainty_plot(uncertainty_mean, title = "Uncertainty mean (BART SDMs)")

# Export uncertainty plot
if (@isdefined save_figures) && save_figures == true
    savefig(plot(unc_plot, dpi = 150), joinpath("fig", "bart", "x_bart_uncertainty.png"))
end
