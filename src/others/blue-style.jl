using JuliaFormatter

# Test formatter on single script
format("./src/04_analysis.jl", BlueStyle())
format("./src/04_analysis.jl", BlueStyle(); pipe_to_function_call=false)

# Use formatter on main scripts
mainfiles = readdir("./src/"; join=true)
filter!(endswith(".jl"), mainfiles)
show(stdout, "text/plain", mainfiles)
problems = []
function formatfiles(files)
    for file in files
        @info "Formatting $(file)"
        try
            format(file, BlueStyle(); pipe_to_function_call=false)
        catch
            @warn "Error while $(file)"
            push!(problems, file)
            continue
        end
        return problems
    end
end
formatfiles(mainfiles)

# Investigate problematic files
format(problems, BlueStyle(); pipe_to_function_call=false)
format(mainfiles[5], BlueStyle(); pipe_to_function_call=false)

# Use formatter on lib scripts
libfiles = readdir("./src/lib/"; join=true)
filter!(endswith(".jl"), libfiles)
formatfiles(libfiles)

# Use formatter on other scripts
otherfiles = readdir("./src/others/"; join=true)
filter!(endswith(".jl"), otherfiles)
formatfiles(otherfiles)
unique(problems)