"""

    verify_jld2_data(path::AbstractString; extract_recent::Bool = false)

Verify if JLD2 files exist and match ZIP archives in `path`. Only the archives are version-controlled
and hosted on the remote repository.
If a JLD2 does not exist, it will be extracted from the archives and written to disk.

A warning will be returned if ZIP archives are more recent than their corresponding JLD2 files,
likely indicating upstream changes. Re-run with `extract_recent` set to `true` to extract the newer files. 
"""
function verify_jld2_data(path::AbstractString; extract_recent::Bool = false)
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
    end
    # Warn if a ZIP archive is more recent than its corresponding JLD2 file
    more_recent = mtime.(zipfiles) .> mtime.(jldfiles)
    if any(more_recent)
        @warn """
              \n One (or more) ZIP archive is more recent than its corresponding JLD2 file.
              \n Run `verify_jld2_data(path; extract_recent = true)` to extract newer version.
              """
        # Extract JLD2 file from more recent archive, if specified
        if extract_recent
            for (i,j) in zip(findall(more_recent), 1:sum(more_recent))
                @info "Reading $(jldfiles[i]) from more recent archive ($(j)/$(sum(more_recent)))"
                _unzip_jld2(zipfiles[i], jldfiles[i])
            end
        end
    end
end

@time verify_jld2_data(joinpath("data", "jld2"))
@time verify_jld2_data(joinpath("data", "jld2"); extract_recent = true)

function _unzip_jld2(zip_path, jld_path)
    zf = ZipFile.Reader(zip_path)
    open(jld_path, "w") do io
        write(io, read(zf.files[1]))
    end
end

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
