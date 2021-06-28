#### 03c - BART predictions ####
include("required.jl")

## Conditional arguments
# save_data = true

## Load predictions from CSV
summaries = CSV.read(joinpath("data", "proc", "bart_summaries.csv"), DataFrame)
varimps = CSV.read(joinpath("data", "proc", "bart_varimps.csv"), DataFrame)
pred_df = CSV.read(
    joinpath("data", "proc", "bart_predictions_prob.csv"), DataFrame; missingstrings=["NA"]
)
lower_df = CSV.read(
    joinpath("data", "proc", "bart_predictions_lower.csv"), DataFrame; missingstrings=["NA"]
)
upper_df = CSV.read(
    joinpath("data", "proc", "bart_predictions_upper.csv"), DataFrame; missingstrings=["NA"]
)
pres_df = CSV.read(
    joinpath("data", "proc", "bart_predictions_pres.csv"), DataFrame; missingstrings=["NA"]
)

## Create Y matrices

# Get matrix Y
Y = replace(Array(pres_df), missing => nothing) |> Array{Union{Nothing,Float32}}
Yprob = replace(Array(pred_df), missing => nothing) |> Array{Union{Nothing,Float32}}
Ylower = replace(Array(lower_df), missing => nothing) |> Array{Union{Nothing,Float32}}
Yupper = replace(Array(upper_df), missing => nothing) |> Array{Union{Nothing,Float32}}
# Set values to nothing if no species present
inds_zeros = _indsnotobs(Y)
Y[inds_zeros, :] .= nothing

## Create distributions

# Load raw distributions (for grid size)
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
spenames = filter(:type => ==("species"), glossary).full_name
distributions = [
    geotiff(SimpleSDMPredictor, joinpath("data", "raster", "distributions_raw.tif"), i) for
    i in eachindex(spenames)
]
raw_distributions = copy(distributions)
# Cut to Quebec coordinates (optional)
# coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
# raw_distributions = [d[coords_qc] for d in raw_distributions]
# Get layer dimensions & limits
dims = size(raw_distributions[1])
lims = boundingbox(raw_distributions[1])

# Create distribution layers
layers = []
for Y in (Y, Yprob, Ylower, Yupper)
    Ydistrib = replace(Y, 0.0 => nothing)
    Ygrids = [Ydistrib[:, col] for col in 1:size(Ydistrib, 2)]
    Ygrids = reshape.(Ygrids, dims...) .|> Array
    local distributions = SimpleSDMResponse.(Ygrids, lims...)
    push!(layers, distributions)
end
distributions, prob_distrib, lower_distrib, upper_distrib = layers;
distributions

## Export results
# save_data = true
if (@isdefined save_data) && save_data == true
    # Export BART distributions
    geotiff(joinpath("data", "raster", "distributions_bart.tif"), distributions)
    # Extras
    geotiff(joinpath("data", "raster", "bart_xtras_prob-distrib.tif"), prob_distrib)
    geotiff(joinpath("data", "raster", "bart_xtras_lower-distrib.tif"), lower_distrib)
    geotiff(joinpath("data", "raster", "bart_xtras_upper-distrib.tif"), upper_distrib)
    # Update placeholder files (as files are too big for version control)
    placeholder_paths = [
        joinpath("data", "raster", "bart_xtras_prob-distrib_placeholder.tif"),
        joinpath("data", "raster", "bart_xtras_lower-distrib_placeholder.tif"),
        joinpath("data", "raster", "bart_xtras_upper-distrib_placeholder.tif"),
    ]
    for p in placeholder_paths
        open(placeholder_path, "w") do io
            write(io, string(Dates.now()))
        end
    end
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])

## Map uncertainty

# Get uncertainty per species
uncertainty = upper_distrib .- lower_distrib
plotSDM2(uncertainty[1]; c=:viridis)

# Uncertainty sum
uncertainty_sum = sum(uncertainty)
plotSDM2(uncertainty_sum; c=:viridis)
# Uncertainty mean
uncertainty_mean = mean(uncertainty)
plotSDM2(uncertainty_mean; c=:viridis)
histogram(uncertainty_mean)

# Plot uncertainty & histogram
function uncertainty_plot(layer; title="")
    unc_map = plotSDM2(
        layer;
        c=:viridis,
        clim=extrema(layer),
        title="Uncertainty",
        colorbar_title="Uncertainty mean",
    )
    unc_hist = histogram(
        layer;
        legend=:none,
        # ylim = maxrange, # xlabel = "Difference",
        title="Distribution of uncertainty values",
        orientation=:horizontal,
    )
    unc_title = plot(; annotation=(0.5, 0.5, "$(title)"), framestyle=:none)
    l = @layout [t{0.01h}; a{0.6w} b{0.38w}]
    unc_plot = plot(unc_title, unc_map, unc_hist; size=(800, 400), layout=l)
    return unc_plot
end
unc_plot = uncertainty_plot(uncertainty_mean; title="Uncertainty mean (BART SDMs)", dpi=200)

# Export uncertainty plot
if (@isdefined save_figures) && save_figures == true
    savefig(unc_plot, joinpath("fig", "bart", "x_bart_uncertainty.png"))
end
