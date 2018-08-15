target_taxon = GBIF.taxon("Cardinalis cardinalis"; rank=:SPECIES)

occ_data = occurrences(target_taxon, Dict{Any,Any}("limit"=>200))
@progress for i in 1:10
    next!(occ_data)
end
showall!(occ_data)

function is_canada_or_us(r::GBIFRecord)
    return r.country ∈ ["Canada", "United States"]
end

qualitycontrol!(occ_data; filters=[have_ok_coordinates, have_both_coordinates, is_canada_or_us])
