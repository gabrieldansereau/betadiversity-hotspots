include("../required.jl")

include("../07_comparison_data.jl")

nspecies = size(raw.Y, 2)
nsites_raw = format(length(raw.richness); commas=true)
nsites_sdm_obs = format(length(sdm.richness); commas=true)
pres_df = CSV.read(
        joinpath("data", "proc", "bart_predictions_pres.csv"), DataFrame; missingstrings=["NA"]
    );
nsites_sdm_tot = nrow(dropmissing(pres_df))

lcbdraw_min, lcbdraw_max = sprintf1.("%.3e", extrema(raw.lcbd))
lcbdsdm_min, lcbdsdm_max = sprintf1.("%.3e", extrema(sdm.lcbd))

bdtot_raw = sprintf1("%.3f", raw.beta_total)
bdtot_sdm = sprintf1("%.3f", sdm.beta_total)

##

R"""
library(SpatialPack)
load(file.path("data", "rdata", "cor-dutilleul.Rdata"))

z_richness_values = list(
    corr = z_richness$corr,
    Fstat = z_richness$Fstat * z_richness$dof, # as reported in summary
    pvalue = z_richness$p.value
)
z_lcbd_values = list(
    corr = z_lcbd$corr,
    Fstat = z_lcbd$Fstat * z_lcbd$dof, # as reported in summary
    pvalue = z_lcbd$p.value
)
"""

@rget z_richness_values z_lcbd_values

sprintf1("%.3f", z_richness_values[:corr])
sprintf1("%.3f", z_richness_values[:Fstat])
sprintf1("%.3e", z_richness_values[:pvalue])
sprintf1("%.3f", z_lcbd_values[:corr])
sprintf1("%.3f", z_lcbd_values[:Fstat])
sprintf1("%.3e", z_lcbd_values[:pvalue])

include("../09_comparison_plots.jl")

richdiff_min, richdiff_max = sprintf1.("%i", extrema(richness_diff))

lcbdcomraw_min, lcbdcomraw_max = sprintf1.("%.3e", extrema(raw.lcbdcom))
lcbdcomsdm_min, lcbdcomsdm_max = sprintf1.("%.3e", extrema(sdm.lcbdcom))

lcbddiff_min, lcbddiff_max = sprintf1.("%.3e", extrema(lcbd_diff))

richres_min, richres_max = sprintf1.("%.3f", extrema(richres_nb_layer))
lcbdres_min, lcbdres_max = sprintf1.("%.3f", extrema(lcbdres_br_layer))


##

include("../05_subareas.jl")

coords_NE = (left=-80.0, right=-60.0, bottom=40.0, top=50.0)
coords_SW = (left=-120.0, right=-100.0, bottom=30.0, top=40.0)

bdtot_NE = sprintf1("%.3f", beta_NE)
bdtot_SW = sprintf1("%.3f", beta_SW)

lcbdne_min, lcbdne_max = sprintf1.("%.3e", extrema(lcbd_NE))
lcbdsw_min, lcbdsw_max = sprintf1.("%.3e", extrema(lcbd_SW))

bdtot_s1, bdtot_s2, bdtot_s3 = sprintf1.("%.3f", beta_values[[1, mid_ind, end]])

lcbdmin_s1 = sprintf1("%.3e", lcbd_mins[1])
lcbdmax_s1 = sprintf1("%.3e", lcbd_maxs[1])
lcbdmin_s2 = sprintf1("%.3e", lcbd_mins[mid_ind])
lcbdmax_s2 = sprintf1("%.3e", lcbd_maxs[mid_ind])
lcbdmin_s3 = sprintf1("%.3e", lcbd_mins[end])
lcbdmax_s3 = sprintf1("%.3e", lcbd_maxs[end])

ascending_threshold_NE = minimum(lcbd_NE)
sprintf1("%.3e", ascending_threshold_NE)
min_indx_NE = findall(x -> x == ascending_threshold_NE, lcbd_NE.grid)
abs_min_NE = median(richness_NE.grid[min_indx_NE])

ascending_threshold_SW = minimum(lcbd_SW)
sprintf1("%.3e", ascending_threshold_SW)
min_indx_SW = findall(x -> x == ascending_threshold_SW, lcbd_SW.grid)
abs_min_SW = median(richness_SW.grid[min_indx_SW])

# medrichmin
# medrichmax
# medlcbdmin
# medlcbdmax
# medgammamin
# medgammamax
# medbdtotmin
# medbttotmax

##

cols = [
    "SCIENTIFIC NAME",
    "LATITUDE",
    "LONGITUDE",
    "OBSERVATION DATE",
    "GROUP IDENTIFIER",
    "APPROVED",
    "SAMPLING EVENT IDENTIFIER"
]

# Load data from CSV files (from file cut with terminal)
@time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), DataFrame; select=cols)

# Prepare data (arrange values & columns)
@time prepare_ebd_data!(df)

# Select subset with specific columns
select!(df, [:species, :latitude, :longitude, :groupIdentifier, :samplingEventIdentifier])

# Remove 1 Aleutian Islands observation with positive longitude
filter!(:longitude => <(0.0), df)

# Remove duplicates (BUT not ok for counts, so-so for dates, see group-observation.jl in src/test/)
n_obs = nrow(df)
df_groups = filter(:groupIdentifier => !ismissing, df)
unique!(df_groups, [:species, :groupIdentifier])
filter!(:groupIdentifier => ismissing, df)
append!(df, df_groups)

n_obs_tot = format(n_obs; commas=true)
n_obs_u = format(nrow(df), commas=true)
n_checklists = format(length(unique(df.samplingEventIdentifier)), commas=true)

# nobs_tot = 27_821_881
# nobs_u = 22_974_330
# nchecklists = 9_103_750
