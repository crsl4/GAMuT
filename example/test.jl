## Used this for the local tests of about 100...feel free to delete

include("example/simulation-one.jl")
for i in 1:100
    run(`julia example/simulation-one.jl $i`)
end 