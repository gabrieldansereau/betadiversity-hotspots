"""
Given a series of values (itr) and a value (p), returns the quantile to which p
is closest. By default, uses `ql=200` steps between 0 and 1.
"""
function find_quantile(itr, p; ql=200)
    qrange = range(0.0; stop=1.0, length=ql)
    q = quantile(filter(.!isnan, itr), qrange)
    return qrange[findmin(abs.(q.-p))[2]]
end

"""
Read a TIFF file and return a matrix
"""
function get_tiff_data(tiff_file)
    # Register GDAL drivers
    GDAL.registerall()

    # Load the dataset
    dataset = GDAL.open(tiff_file, GDAL.GA_ReadOnly)

    # Band
    band = GDAL.getrasterband(dataset, 1)

    # Matrix
    xs = GDAL.getrasterxsize(dataset)
    ys = GDAL.getrasterysize(dataset)

    bandtype = GDAL.getrasterdatatype(band)

    V = zeros(Float64, (xs, ys))

    GDAL.rasterio(
        band,
        GDAL.GF_Read,
        0, 0, xs, ys,
        pointer(V),
        xs, ys,
        GDAL.getrasterdatatype(band),
        0, 0
        )

    K = zeros(Float64, (ys, xs))
    for (i,r) in enumerate(reverse(1:size(V, 2)))
        K[i,:] = V[:,r]
    end

    this_min = minimum(V)

    for i in eachindex(K)
        K[i] = K[i] > this_min ? K[i] : NaN
    end

    return K

end

function get_closest_grid_point(x::NTuple{2,Float64}, lon::Vector{Float64}, lat::Vector{Float64})
    dx = findmin(abs.(x[1].-lon))[2]
    dy = findmin(abs.(x[2].-lat))[2]
    return (dx, dy)
end
