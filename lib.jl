"""
Given a series of values (itr) and a value (p), returns the quantile to which p
is closest. By default, uses `ql=200` steps between 0 and 1.
"""
function find_quantile(itr, p; ql=100)
    qrange = range(0.0; stop=1.0, length=ql)
    q = quantile(itr, qrange; sorted=true)
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
    return (dy, dx)
end

function get_value_at_position(x, v, lon, lat)
    tx, ty = get_closest_grid_point(x, lon, lat)
    return v[tx, ty]
end


function get_bioclim_values(records::GBIFRecords, bioclim_variables, longitudes, latitudes)
    n_variables = length(bioclim_variables)
    n_observations = sum(records.show)
    values_table = Vector{Vector{Float64}}()
    for (v_index, bioclim_variable) in enumerate(bioclim_variables)
        temp = [get_value_at_position((record.longitude, record.latitude), bioclim_variable, longitudes, latitudes) for record in records]
        these_values = sort(filter(.!isnan, temp))
        push!(values_table, these_values)
    end
    return values_table
end


function get_quantile_matrix(bioclim_variable, observations, b_box, longitudes, latitudes)
    b_pts = [get_closest_grid_point(b, longitudes, latitudes) for b in b_box]
    x_span = b_pts[1][1]:b_pts[2][1]
    y_span = b_pts[1][2]:b_pts[2][2]

    q_matrix = zeros(Float64, (length(x_span), length(y_span)))

    for (i, x) in enumerate(x_span)
        for (j, y) in enumerate(y_span)
            local_val = bioclim_variable[x, y]
            if isnan(local_val)
                q_matrix[i,j] = NaN
            else
                this_q = find_quantile(observations, local_val)
                this_q = this_q > 0.5 ? 1.0-this_q : this_q
                q_matrix[i,j] = this_q
            end

        end
    end
    return 2.0 .* q_matrix
end
