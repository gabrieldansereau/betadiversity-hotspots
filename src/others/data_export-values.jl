import Pkg; Pkg.activate(".")
include("../required.jl")

include("../07a_comparison_data.jl")

nspecies = size(raw.Y, 2)
nsites_raw = length(raw.richness)
nsites_sdm = length(sdm.richness)
coords_full = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
nsites_tot = length(worldclim(1)[coords_full])

lcbdraw_min, lcbdraw_max = extrema(raw.lcbd)
lcbdsdm_min, lcbdsdm_max = extrema(sdm.lcbd)

lcbdcomraw_min, lcbdcomraw_max = extrema(raw.lcbdcom)
lcbdcomsdm_min, lcbdcomsdm_max = maximum(sdm.lcbdcom)

bdtot_raw = raw.beta_total
bdtot_sdm = sdm.beta_total

##

include("../07c_comparison_plots.jl")

richdiff_min, richdiff_max = extrema(richness_diff)
lcbddiff_min, lcbddiff_max = extrema(lcbd_diff)

richres_min, richres_max = extrema(richres_nb_layer)
lcbdres_min, lcbdres_max = extrema(lcbdres_br_layer)

##

include("../05_subareas.jl")

bdtot_NE = beta_NE
bdtot_SW = beta_SW

lcbdne_min, lcbdne_max = extrema(lcbd_NE)
lcbdsw_min, lcbdsw_max = extrema(lcbd_SW)

bdtot_s1, bdtot_s2, bdtot_s3 = beta_values[[1, mid_ind, end]]

lcbdmin_s1 = lcbd_mins[1]
lcbdmax_s1 = lcbd_maxs[1]
lcbdmin_s2 = lcbd_mins[2]
lcbdmax_s2 = lcbd_maxs[2]
lcbdmin_s3 = lcbd_mins[3]
lcbdmax_s3 = lcbd_maxs[3]

# medrichmin
# medrichmax
# medlcbdmin
# medlcbdmax
# medgammamin
# medgammamax
# medbdtotmin
# medbttotmax

##

@time df = CSV.read(joinpath("data", "proc", "ebd_warblers_prep.csv"), DataFrame, header=true, delim="\t");

nobs_tot = 27_821_881
nobs_u = 22_974_330
nchecklists =  9_103_750

