function worldclim(i::Int64; path="assets")
    @assert 1 ≤ i ≤ 19
    code = lpad(i, 2, "0")
    raw_data = get_data_from_tiff(joinpath(path, "wc2.0_bio_10m_$(i).tif"))
    return SDMPredictor(raw_data, -180.0, 180.0, -90.0, 90.0)
end
