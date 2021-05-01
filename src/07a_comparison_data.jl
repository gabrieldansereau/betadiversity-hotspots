import Pkg; Pkg.activate(".")
include("required.jl")

## Run analysis on raw data
outcome = "raw"
include("04_analysis.jl")

# Add non-relative LCBD values
lcbdnr = calculate_lcbd(Y, lcbd; relative = false)

# Assemble results
raw = (Y = Y,
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
sdm = (Y = Y,
       richness = richness,
       richness_plot = richness_plot,
       lcbd = lcbd,
       lcbdnr = lcbdnr,
       beta_total = beta_total,
       lcbd_plot = lcbdtr_plot,
       rel_plot = rel2d_plot
       )

# Export to JLD2
jld_path = joinpath("data", "jld2", "comparison-results.jld2")
@save jld_path raw sdm

## Prepare data for GLMs
# Assemble in DataFrame
results = DataFrame([raw.richness, raw.lcbd, sdm.richness, sdm.lcbd, raw.lcbdnr, sdm.lcbdnr])
rename!(results, 
        :x1 => :richness_raw, :x2 => :lcbd_raw, 
        :x3 => :richness_sdm, :x4 => :lcbd_sdm,
        :x5 => :lcbdnr_raw, :x6 => :lcbdnr_sdm)
# Remove lines with missing values
allowmissing!(results, Not([:latitude, :longitude]))
for col in eachcol(results)
    replace!(col, nothing => missing)
end
dropmissing!(results)
# Convert to Float64 (cannot write to CSV otherwise)
results = convert.(Float64, results)

# Export to CSV
CSV.write(joinpath("data", "proc", "comparison-results.csv"), results, delim = "\t")
