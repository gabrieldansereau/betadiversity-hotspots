function worldclim(i::Int64; resolution="5", path="assets")
    @assert 1 ≤ i ≤ 19
    code = lpad(i, 2, "0")
    raw_data = get_data_from_tiff(joinpath(path, "wc2.0_bio_$(resolution)m_$(code).tif"))
    return SimpleSDMPredictor(raw_data, -180.0, 180.0, -90.0, 90.0)
end
