using Plots
using Statistics
using Distributions

gridsize = 100
x = range(0.0; stop=1.0, length=gridsize)
y = range(0.0; stop=1.0, length=gridsize)
S = zeros(Float64, (gridsize, gridsize))

function find_quantile(itr, p)
    qrange = range(0.0; stop=1.0, length=200)
    q = quantile(itr, qrange)
    return qrange[findmin(abs.(q.-p))[2]]
end

Dx = TruncatedNormal(0.4, 0.2, 0.1, 0.6)
Dy = TruncatedNormal(0.6, 0.1, 0.2, 0.8)

nx = rand(Dx, 105)
ny = rand(Dy, 105)

for (i, ax) in enumerate(x)
    qx = find_quantile(nx, ax)
    qx = qx > 0.5 ? 1-qx : qx
    for (j, ay) in enumerate(y)
        qy = find_quantile(ny, ay)
        qy = qy > 0.5 ? 1-qy : qy
        this_q = minimum([qx, qy])
        S[i,j] = this_q
        S[i,j] = this_q > 0.5 ? 1 - this_q : this_q
    end
end

S = S .* 2.0

heatmap(x, y, S, c=:Blues, aspectratio=1, frame=:tick)
scatter!(ny, nx, lab="", msw=0)
yaxis!((0,1))
xaxis!((0,1))
#
