function _unzip_jld2(zip_path, jld_path)
    zf = ZipFile.Reader(zip_path)
    open(jld_path, "w") do io
        write(io, read(zf.files[1]))
    end
end
function verify_jld2_data(path)
    # List zip files
    files = readdir(path, join = true)
    zipfiles = filter(x -> occursin(".zip", x), files)
    # Check if corresponding .jld2 exist
    jldfiles = replace.(zipfiles, ".zip" => ".jld2")
    missing_files = map(!isfile, jldfiles)
    # Extract missing files from archive & write to disk
    if any(missing_files)
        for i in findall(missing_files)
            @info "Reading $(jldfiles[i]) from archive ($(i)/$(sum(missing_files)))"
            _unzip_jld2(zipfiles[i], jldfiles[i])
        end
    else
        @info "All JLD2 files exist"
    end
end
@time verify_jld2_data(joinpath("data", "jld2"))

function _zip_jld2(zip_path, jld_path)
    w = ZipFile.Writer(zip_path) # create archive folder
    filename = split(jld_path, "/")[end]
    f = ZipFile.addfile(w, filename, method=ZipFile.Deflate) # create file for compression
    open(jld_path, "r") do io # file to compress
        write(f, io) # compress file, ~60sec
    end
    close(w)
end
@time _zip_jld2(zipfiles[2], jldfiles[2])
