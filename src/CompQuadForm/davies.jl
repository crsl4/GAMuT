function davies(q::Float64, lambda::Vector) 
    h = fill(1, length(lambda))
    df = convert.(Int64, h)
    delta = zeros(length(lambda))
    noncentral = convert.(Float64, delta)
    sigma = 0.0
    lim = 10000
    acc = 0.0001
    r = convert(Int64, length(lambda))
    s = r
    trace = fill(1, 7) #was zeros in source but order function needs 1's for indexing 
    ifault = 0
    res = 0.0

    if any(delta .< 0)
        error("All non centrality parameters in 'delta' should be positive!")
    end
    if length(h) != r
        error("lambda and h should have the same length!")
    end
    if length(delta) != r
        error("lambda and delta should have the same length!")
    end
    include("bin/GAMuT/src/CompQuadForm/c_davies.jl")
    out = qfc(vec(lambda), noncentral, df, s, r, sigma, q, lim, acc, ifault, res, 0.0, vec(trace), true)
    out[:res] = 1 - out[:res]
    if out[:res] > 1 
        warn("Consider playing with 'lim' or 'acc'.")
    end
    return([trace=out[:trace], ifault=out[:ifault], Qq=out[:res]])
end
