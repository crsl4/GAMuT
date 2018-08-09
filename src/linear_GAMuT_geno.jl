#Genotypic Similarity: linear kernel
function linear_GAMuT_geno(x)
    lc = x'*x #Made specifically for eigvals function
    #https://stackoverflow.com/questions/51731339/eigenvalues-of-a-matrix-assuming-symmetry
    ev_Lc = eigvals(Symmetric(lc, :L))
    lc = A_mul_Bt(x, x)
    ev_Lc = ev_Lc[ev_Lc .> 1e-08]
    sort!(ev_Lc,rev=true)
    return([lc,ev_Lc])
end
