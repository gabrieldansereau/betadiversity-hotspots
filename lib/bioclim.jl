function bioclim(p::SDMLayer, r::GBIFRecords)
    observed_values = p[r]
    qfinder = ecdf(observed_values)
    lq = zeros(Float64, size(p))
    for i in eachindex(p.grid)
        if isnan(p.grid[i])
            local_quantile = NaN
        else
            local_quantile = qfinder(p.grid[i])
            if local_quantile > 0.5
                local_quantile = 1.0-local_quantile
            end
            local_quantile = 2.0 * local_quantile
        end
        lq[i] = local_quantile
    end
    prediction = SDMLayer(lq, p.left, p.right, p.bottom, p.top)
end
