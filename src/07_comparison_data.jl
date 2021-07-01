include("required.jl")

## Conditional argument
# save_additional_data = false
# save_additional_data = true

## Run analysis on raw data
outcome = "raw"
include("04_full-extent.jl")

# Assemble results
raw = (
    Y=Y,
    richness=richness_full,
    richness_plot=richness_plot,
    lcbd=lcbd_full,
    beta_total=beta_full,
    lcbd_plot=lcbd_plot,
    rel_plot=rel2d_plot,
)

## Run analysis on sdm data
outcome = "bart"
include("04_full-extent.jl")

# Assemble results
sdm = (
    Y=Y,
    richness=richness_full,
    richness_plot=richness_plot,
    lcbd=lcbd_full,
    beta_total=beta_full,
    lcbd_plot=lcbd_plot,
    rel_plot=rel2d_plot,
)

## Recalculate LCBD values on same sites

# Get indices in common (with values for both raw & sdm)
inds_common = intersect(_indsobs(raw.Y), _indsobs(sdm.Y))

# Recalculate LCBD values on same sites
function common_lcbd(Y, layer, inds_common)
    Ycommon = copy(Y)
    Ycommon[Not(inds_common), :] .= nothing
    lcbd_common = lcbd(Ycommon, layer; relative=false)
    return lcbd_common
end
lcbdcom_raw = common_lcbd(raw.Y, raw.lcbd, inds_common)
lcbdcom_sdm = common_lcbd(sdm.Y, sdm.lcbd, inds_common)

# Add to NamedTuples
raw = (raw..., lcbdcom=lcbdcom_raw)
sdm = (sdm..., lcbdcom=lcbdcom_sdm)

# Export to JLD2
if (@isdefined save_additional_data) && save_additional_data == true
    # Export to JLD2
    jld_path = joinpath("data", "jld2", "comparison-results.jld2")
    @save jld_path raw sdm
    # Update placeholder file (as file is too big for version control)
    placeholder_path = joinpath("data", "jld2", "comparison-results_placeholder.jld2")
    open(placeholder_path, "w") do io
        write(io, string(Dates.now()))
    end
end

## Prepare data for GLMs

# Assemble in DataFrame
results = DataFrame([raw.richness, raw.lcbdcom, sdm.richness, sdm.lcbdcom])
rename!(
    results, :x1 => :richness_raw, :x2 => :lcbd_raw, :x3 => :richness_sdm, :x4 => :lcbd_sdm,
)

# Remove lines with missing values
dropmissing!(results)

# Convert to Float64 (cannot write to CSV otherwise)
results = convert.(Float64, results)

# Export to CSV
if (@isdefined save_additional_data) && save_additional_data == true
    CSV.write(joinpath("data", "proc", "comparison-results.csv"), results; delim="\t")
end
