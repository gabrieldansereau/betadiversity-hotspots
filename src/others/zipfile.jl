import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))
using ZipFile
cd("../test/data/")

# Load original file
@load "rf-distributions.jld2" distributions
rfdist = distributions

# Unzip existing zip file
zf = ZipFile.Reader("rf-distributions.zip")
@time write("test.jld2", read(zf.files[1]))
# Load result
@load "test.jld2" distributions
zfdist = copy(distributions)

# Zip existing file
w = ZipFile.Writer("example.zip")
f = ZipFile.addfile(w, "test2.jld2", method=ZipFile.Deflate)
@time write(f, open("test.jld2", "r"))
close(w)
# Unzip result & load
zf2 = ZipFile.Reader("example.zip")
@time write("test2.jld2", read(zf2.files[1]))
@load "test2.jld2" distributions
zfdist2 = copy(distributions)

# Compare results
map(x -> filter(!isnan, x.grid), zfdist) == map(x -> filter(!isnan, x.grid), rfdist)
map(x -> filter(!isnan, x.grid), zfdist2) == map(x -> filter(!isnan, x.grid), rfdist)

## Zipfile documentation
using ZipFile
w = ZipFile.Writer("example.zip");
f = ZipFile.addfile(w, "hello.txt");
write(f, "hello world!")
f = ZipFile.addfile(w, "julia.txt", method=ZipFile.Deflate);
write(f, "Julia");
close(w)
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
