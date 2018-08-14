include("simulation_one_files/simulation-one.jl")
for i in 1:1000
    run(`julia simulation_one_files/simulation-one.jl $i`)
end 