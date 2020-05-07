import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))
using ZipFile
cd("../test/data/")

# Load original file
@load "rf-distributions.jld2" distributions
rfdist = distributions

## Zip existing file
@time begin
    w = ZipFile.Writer("$(outcome)-distributions.zip") # create archive folder
    f = ZipFile.addfile(w, "$(outcome)-distributions.jld2", method=ZipFile.Deflate) # create file for compression
    open("$(outcome)-distributions.jld2", "r") do io # file to compress
        write(f, io) # compress file, ~60sec
    end
    close(w)
end

## Unzip a file

# Unzip previously archived file & write to disk
zf = ZipFile.Reader("$(outcome)-distributions.zip")
@time write("unziptest.jld2", read(zf.files[1]))
@load "unziptest.jld2" distributions
zfdist = copy(distributions)

# Compare results
map(x -> filter(!isnan, x.grid), zfdist) == map(x -> filter(!isnan, x.grid), rfdist)

## Zipfile documentation
using ZipFile
# Zip files
w = ZipFile.Writer("example.zip");
f = ZipFile.addfile(w, "hello.txt");
write(f, "hello world!")
f = ZipFile.addfile(w, "julia.txt", method=ZipFile.Deflate);
write(f, "Julia");
close(w)
# Unzip files
r = ZipFile.Reader("example.zip");
for f in r.files
      println("Filename: $(f.name)")
      write(stdout, read(f, String));
end
close(r)

## Zip extract from SimpleSDMLayers
for rf in zf.files
    if joinpath(path, rf.name) in paths
        if !isfile(joinpath(path, rf.name))
            @info "Reading layer $(rf.name) from archive"
            write(joinpath(path, rf.name), read(rf))
        end
    end
end
close(zf)
