## Shapefiles functions useful for mapping

# Download worldmap shapefile
function download_shapefile(res)
    @assert res ∈ [50, 100]
    dir =
        "https://github.com/nvkelso/natural-earth-vector/" * "raw/master/$(res)m_physical/"

    fn = "ne_$(res)m_land.shp"
    return run(`wget $dir/$fn -P ./assets/`)
end

# Read worlmap shapefile
function worldshape(res)
    handle = open("./assets/ne_$(res)m_land.shp", "r") do io
        read(io, Shapefile.Handle)
    end
    return handle
end

# Clip worldmap shapefile
function clip(s::Shapefile.Handle, l::SimpleSDMLayer)
    return filter(x -> isin(x, l), s.shapes)
end

# Check layer bounds are outside polygon
function isin(p::Shapefile.Polygon, l::SimpleSDMLayer)
    out = false
    for xy in p.points
        xy.x < l.right && (out = true)
        xy.x > l.left && (out = true)
        xy.y < l.top && (out = true)
        xy.y > l.bottom && (out = true)
    end
    return out
end
