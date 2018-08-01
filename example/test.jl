include("gamut_simulation/simulation-one.jl")
for i in 1:1000
    run(`julia gamut_simulation/simulation-one.jl $i`)
end 