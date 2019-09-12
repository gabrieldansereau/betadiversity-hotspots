using Distributed
using JLD2
using RCall

@time @everywhere include("src/required.jl")

df = CSV.read("../data/ebd/processed/ebd_warblers_cut.csv", header=true, delim="\t")
head = CSV.read("../data/ebd/ebd_warblers_head.csv", header=true, delim="\t")
prep = CSV.read("../data/ebd/ebd_warblers_prep.csv", header=true, delim="\t")
@load "../data/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred
@load "../data/pres-abs-ebd.jld2" pres_abs
wc_vars = pmap(worldclim, 1:19)

# Prepare data (arrange values & columns)
df = prepare_ebd_data(df)
head = prepare_ebd_data(head)

# Remove duplicates (BUT not ok for counts, so-so for dates, see group-observation.jl in src/test/)
df_nogroups = filter(x -> ismissing(x[:groupIdentifier]), df)
df_groups = dropmissing(df, :groupIdentifier)
df_groups_unique = unique(df_groups, [:species, :groupIdentifier])
newdf = vcat(df_nogroups, df_groups_unique)

head_nogroups = filter(x -> ismissing(x[:groupIdentifier]), head)
head_groups = dropmissing(head, :groupIdentifier)
head_groups_unique = unique(head_groups, [:species, :groupIdentifier])
newhead = vcat(head_nogroups, head_groups_unique)

$library(xtable)

# Variables base
varbase = select(newhead, [:species, :latitude, :longitude, :observationDate, :allSpeciesReported, :observationCount])
@rput varbase
$varbase$observationDate = as.character(varbase$observationDate)
$xtable(head(varbase), format="latex")

# Variable échant
varechant = select(newhead, [:protocolType, :durationMinutes, :effortDistanceKm, :numberObservers])
@rput varechant
$xtable(head(varechant), format="latex")

speciescounts = by(newdf, :species, nrow)
protocolcounts = by(newdf, :protocolType, nrow)
show(sort(protocolcounts, :x1, rev=true), allrows=true)
