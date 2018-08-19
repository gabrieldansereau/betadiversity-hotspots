using Shapefile

for res in [50,110]
    dir = "https://github.com/nvkelso/natural-earth-vector/" *
        "raw/master/$(res)m_physical/"

    fn = "ne_$(res)m_land.shp"
    run(`wget $dir/$fn -P ./assets/`)
end
