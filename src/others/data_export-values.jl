include("../required.jl")

include("../07a_comparison_data.jl")

nspecies = size(raw.Y, 2)
nsites_raw = format(length(raw.richness); commas=true)
nsites_sdm = format(length(sdm.richness); commas=true)
coords_full = (left=-145.0, right=-50.0, bottom=20.0, top=75.0)
nsites_tot = format(
    length(SimpleSDMPredictor(WorldClim, BioClim, 1; coords_full...)); commas=true
)

lcbdraw_min, lcbdraw_max = sprintf1.("%.3e", extrema(raw.lcbd))
lcbdsdm_min, lcbdsdm_max = sprintf1.("%.3e", extrema(sdm.lcbd))

lcbdcomraw_min, lcbdcomraw_max = sprintf1.("%.3e", extrema(raw.lcbdcom))
lcbdcomsdm_min, lcbdcomsdm_max = sprintf1.("%.3e", extrema(sdm.lcbdcom))

bdtot_raw = sprintf1("%.3f", raw.beta_total)
bdtot_sdm = sprintf1("%.3f", sdm.beta_total)

##

include("../07c_comparison_plots.jl")

richdiff_min, richdiff_max = sprintf1.("%i", extrema(richness_diff))
lcbddiff_min, lcbddiff_max = sprintf1.("%.3e", extrema(lcbd_diff))

richres_min, richres_max = sprintf1.("%.3f", extrema(richres_nb_layer))
lcbdres_min, lcbdres_max = sprintf1.("%.3f", extrema(lcbdres_br_layer))

R"""
library(SpatialPack)
load(file.path("data", "rdata", "cor-dutilleul.Rdata"))
z_richness
z_richness$p.value
z_lcbd
z_lcbd$p.value
"""
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

@time df = CSV.read(
    joinpath("data", "proc", "ebd_warblers_prep.csv"), DataFrame; header=true, delim="\t"
);

nobs_tot = 27_821_881
nobs_u = 22_974_330
nchecklists = 9_103_750
