"""
Given a series of values (itr) and a value (p), returns the quantile to which p
is closest. By default, uses `ql=200` steps between 0 and 1.
"""
function find_quantile(itr, p; ql=200)
    qrange = range(0.0; stop=1.0, length=ql)
    q = quantile(itr, qrange)
    return qrange[findmin(abs.(q.-p))[2]]
end
