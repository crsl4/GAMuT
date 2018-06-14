using RCall
R"""
source("all_gamut_functions.R")
source("simulation-one.r")


b = replicate(1, rnorm(5000))


x = linear_GAMuT_geno(b)
y = TestGAMuT(x$Lc,x$ev_Lc,x$Lc,x$ev_Lc)
"""
@rget y 
@rget b
@rget a
@rget c

include("linear_GAMuT_geno.jl")
include("proj_GAMuT_pheno.jl")
include("test_GAMuT.jl")
lc,ev_Lc = linear_GAMuT_geno(b)
y2 = testGAMuT(lc,ev_Lc,lc,ev_Lc)

isapprox(y,y2) || error("julia and R functions do not match")

R"""
z = proj_GAMuT_pheno(b)
w = TestGAMuT(z$Kc, z$ev_Kc, z$Kc, z$ev_Kc)
"""

@rget w
kc,ev_Kc = proj_GAMuT_pheno(b)
w2 = testGAMuT(kc,ev_Kc,kc,ev_Kc)

isapprox(w, w2) || error("julia and R functions do not match")