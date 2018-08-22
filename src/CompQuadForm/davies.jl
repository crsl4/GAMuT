function davies(q::Float64, lambda::Vector) 
    h = fill(1, length(lambda))
    df = convert.(Int64, h)
    delta = zeros(length(lambda))
    noncentral = convert.(Float64, delta)
    sigma = 0.0
    tausq = 0.0
    lim = 10000
    acc = 0.0001
    r = convert(Int64, length(lambda))
    s = r
    trace = fill(1, 8) #for indexing (C is 0 based so goes from 0-7)
    ifault = 0
    res = 0.0
    th = vec(lambda)
    ndstart = true

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
    
    ## th is not listed in source (it is an external vector), and I do not know what the values in th are supposed to be so I set it equal to lambda 
    ## also do not know what ndstart = so I set it to true so it will go through order function 
    out = qfc(vec(lambda), noncentral, df, r, sigma, q, lim, acc, trace, ifault, res, s, tausq, th, ndstart)
    out[:res] = 1 - out[:res]
    if out[:res] > 1 
        warn("Consider playing with 'lim' or 'acc'.")
    end
    return([trace=out[:trace], ifault=out[:ifault], Qq=out[:res]])
end
