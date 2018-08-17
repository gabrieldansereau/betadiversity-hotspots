function find_quantile(itr, p; ql=150)
    qrange = range(0.0; stop=1.0, length=ql)
    q = quantile(itr, qrange; sorted=true)
    return qrange[findmin(abs.(q.-p))[2]]
end

function bioclim(p::SDMLayer, r::GBIFRecords)
    observed_values = sort(p[r])
    lq = zeros(Float64, size(p))
    for i in eachindex(p.grid)
        if isnan(p.grid[i])
            local_quantile = NaN
        else
            local_quantile = find_quantile(observed_values, p.grid[i])
            if local_quantile > 0.5
                local_quantile = 1.0-local_quantile
            end
            local_quantile = 2.0 * local_quantile
        end
        lq[i] = local_quantile
    end
    prediction = SDMLayer(lq, p.left, p.right, p.bottom, p.top)
end
