using Plots

x = sort(rand(1:10, 10))
y = sort(rand(1:20, 10))

plot(x,y)

plot(x, y; title = "test", color = "black")

function plotthis(x, y; kw...)
    plot(x, y; kw...)
end # this way works

plotthis(x, y)
plotthis(x, y; title = "test")
plotthis(x, y; title = "test", color = "black")

function plotthat(x, y; args...)
    plot(x, y; args...)
end # keyword varargs can a different name

plotthat(y, x, color = "black") # ; not needed, as always for keyword arguments
plotthat(y, x; :color => "black") # symbols also work
plotthat(y, x, :color => "red", :style => :dash) # and symbols don't require semi-colons either

function plotit(x, y, args...)
    plot(x, y, args...)
end

plotit(x,y) # works
plotit(x,y, color = "black") # doesn't accept keyword arguments
plotit(x,y, "black") # not interpreted as color

function plotwhat(x, y; title)
    plot(x, y; title = title)
end

plotwhat(x, y) # keyword arguments have to be assigned

## Summary

#=
varargs: f(x, y, ...)
optional: f(x, y, z = 1)
keywords: f(x, y; color) # semi-colon required at definition



=#
