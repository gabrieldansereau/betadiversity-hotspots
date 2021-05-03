import Pkg; Pkg.activate(".")
include("required.jl")

## Run analysis on raw data
outcome = "raw"
include("04_analysis.jl")

# Add non-relative LCBD values
lcbdnr = calculate_lcbd(Y, lcbd; relative = false)

# Assemble results
raw = (
    Y = Y,
    richness = richness,
    richness_plot = richness_plot,
    lcbd = lcbd,
    lcbdnr = lcbdnr,
    beta_total = beta_total,
    lcbd_plot = lcbdtr_plot,
    rel_plot = rel2d_plot
)

## Run analysis on sdm data
outcome = "bart"
include("04_analysis.jl")

# Add non-relative LCBD values
lcbdnr = calculate_lcbd(Y, lcbd; relative = false)

# Assemble results
sdm = (
    Y = Y,
    richness = richness,
    richness_plot = richness_plot,
    lcbd = lcbd,
    lcbdnr = lcbdnr,
    beta_total = beta_total,
    lcbd_plot = lcbdtr_plot,
    rel_plot = rel2d_plot
)

## Recalculate LCBD values on same sites
inds_common = intersect(_indsobs(raw.Y), _indsobs(sdm.Y))
function recalculate_lcbd(Y, layer, inds_common)
    Ycommon = copy(Y)
    Ycommon[Not(inds_common), :] .= nothing
    lcbd_common = calculate_lcbd(Ycommon, layer; relative = false)
    return lcbd_common
end
lcbdcom_raw = recalculate_lcbd(raw.Y, raw.lcbd, inds_common)
lcbdcom_sdm = recalculate_lcbd(sdm.Y, sdm.lcbd, inds_common)
plotSDM2(lcbdcom_raw, c = :viridis)
plotSDM2(lcbdcom_sdm, c = :viridis)
extrema(lcbdcom_raw)
extrema(lcbdcom_sdm)
extrema(raw.lcbd)
extrema(sdm.lcbd)

# Add to NamedTuples
raw = (raw..., lcbdcom = lcbdcom_raw)
sdm = (sdm..., lcbdcom = lcbdcom_sdm)
extrema(raw.lcbdcom)
extrema(sdm.lcbdcom)
extrema(raw.lcbd)
extrema(sdm.lcbd)

## Tests
# isequal([sdm.richness[i] for i in inds_common], results.richness_sdm)
# isequal([raw.richness[i] for i in inds_common], results.richness_raw)
# isequal([sdm.lcbd[i] for i in inds_common], results.lcbd_sdm)
# isequal([raw.lcbd[i] for i in inds_common], results.lcbd_raw)
# Ycommon = copy(raw.Y)
# Ycommon[Not(inds_common), :] .= nothing
# isequal(collect(calculate_richness(Ycommon, raw.richness)), raw.richness.grid[inds_common])
# calculate_lcbd(Ycommon, raw.lcbd; relative = false) |>
#     x -> plotSDM2(x, c = :viridis)
# plotSDM2(raw.lcbd; relative = false, c = :viridis)

# Export to JLD2
jld_path = joinpath("data", "jld2", "comparison-results.jld2")
@save jld_path raw sdm

## Prepare data for GLMs
# Assemble in DataFrame
results = DataFrame([raw.richness, raw.lcbdcom, sdm.richness, sdm.lcbdcom, raw.lcbdnr, sdm.lcbdnr])
rename!(
    results, 
    :x1 => :richness_raw, :x2 => :lcbd_raw, 
    :x3 => :richness_sdm, :x4 => :lcbd_sdm,
    :x5 => :lcbdnr_raw, :x6 => :lcbdnr_sdm
)
# Remove lines with missing values
allowmissing!(results, Not([:latitude, :longitude]))
for col in eachcol(results)
    replace!(col, nothing => missing)
end
dropmissing!(results)
# Convert to Float64 (cannot write to CSV otherwise)
results = convert.(Float64, results)
extrema(results.lcbd_raw)
extrema(results.lcbd_sdm)
extrema(results.lcbdnr_raw)
extrema(results.lcbdnr_sdm)

# Export to CSV
CSV.write(joinpath("data", "proc", "comparison-results.csv"), results, delim = "\t")
