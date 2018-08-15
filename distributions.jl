target_taxon = GBIF.taxon("Haemorhous purpureus"; rank=:SPECIES)

finch_occ = occurrences(target_taxon)
@progress for i in 1:50
    next!(finch_occ)
end
showall!(finch_occ)

function is_canada_or_us(r::GBIFRecord)
    return r.country ∈ ["Canada", "United States"]
end

qualitycontrol!(finch_occ; filters=[have_ok_coordinates, have_both_coordinates, is_canada_or_us])

lat = []
lon = []

for i in finch_occ
    @info (latitude = i.latitude, longitude = i.longitude)
    push!(lat, i.latitude)
    push!(lon, i.longitude)
end

writedlm("assets/$(replace(target_taxon.name, " " => "_")).coordinates", hcat(lon, lat))
