using RCall
R"""
source("all_gamut_functions.R")
"""
a = [1 2 3; 4 5 6]
b = [4 9; 5 10; 8 11]
c = [1 4; 2 5; 3 6]
d = [0.182378 0.764705 .351618; .599969 .786265 .0844763; .43153 .969674 .408121; .233296 .250268 .641506]
R"""
d = matrix(c(0.182378, 0.764705, 0.351618, .599969, .786265, .0844763, .43153, .969674, .408121, .233296, .250268, .641506), nrow=4, ncol=3, byrow=TRUE)
x = linear_GAMuT_geno(d)
y = TestGAMuT(x$Lc,x$ev_Lc,x$Lc,x$ev_Lc)
"""
@rget y 

include("linear_GAMuT_geno.jl")
include("proj_GAMuT_pheno.jl")
include("test_GAMuT.jl")
lc,ev_Lc = linear_GAMuT_geno(d)
y2 = testGAMuT(lc,ev_Lc,lc,ev_Lc)

isapprox(y,y2) || error("julia and R functions do not match")

R"""
z = proj_GAMuT_pheno(d)
w = TestGAMuT(z$Kc, z$ev_Kc, z$Kc, z$ev_Kc)
"""

@rget w
kc,ev_Kc = proj_GAMuT_pheno(d)
w2 = testGAMuT(kc,ev_Kc,kc,ev_Kc)

isapprox(w, w2) || error("julia and R functions do not match")
