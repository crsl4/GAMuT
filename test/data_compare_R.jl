using RCall
R"""
source("src/all_gamut_functions.r")
source("test/simulation-one-edited.r")

a = UNR_OBS

x = linear_GAMuT_geno(a)
y = TestGAMuT(x$Lc,x$ev_Lc,x$Lc,x$ev_Lc)
"""
@rget y 
@rget x
@rget a

include("src/linear_GAMuT_geno.jl")
include("src/proj_GAMuT_pheno.jl")
include("src/test_GAMuT.jl")
lc,ev_Lc = linear_GAMuT_geno(a)
y2 = testGAMuT(lc,ev_Lc,lc,ev_Lc)

all(isapprox.(x[:Lc],lc)) || error("julia and R functions do not match")
all(isapprox.(x[:ev_Lc],ev_Lc)) || error("julia and R functions do not match")
isapprox(y,y2) || error("julia and R functions do not match")


R"""
##Need square matrix
b = replicate(500, rnorm(500))

z = proj_GAMuT_pheno(b)
w = TestGAMuT(z$Kc, z$ev_Kc, z$Kc, z$ev_Kc)
"""
@rget b
@rget w

kc,ev_Kc = proj_GAMuT_pheno(b)
w2 = testGAMuT(kc,ev_Kc,kc,ev_Kc)

isapprox(w, w2) || error("julia and R functions do not match")