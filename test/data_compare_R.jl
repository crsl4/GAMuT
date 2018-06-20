using RCall
##------------------------------------------------------------------------------
## Linear Kernel: 
##------------------------------------------------------------------------------
R"""
source("./src/all_gamut_functions.r")
source("./src/simulation-one.r")

x1 = linear_GAMuT_geno(Y1)
x2 = linear_GAMuT_geno(Y2)
y = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""
@rget Y1
@rget Y2
@rget x1
@rget x2
@rget y 

include("./src/linear_GAMuT_geno.jl")
include("./src/proj_GAMuT_pheno.jl")
include("./src/test_GAMuT.jl")

lc,ev_Lc = linear_GAMuT_geno(Y1)
lc_2,ev_Lc_2 = linear_GAMuT_geno(Y2)
y2 = testGAMuT(lc,ev_Lc,lc_2,ev_Lc_2)

all(isapprox.(x1[:Lc],lc)) || error("julia and R functions do not match")
all(isapprox.(x1[:ev_Lc],ev_Lc)) || error("julia and R functions do not match")
all(isapprox.(x2[:Lc],lc_2)) || error("julia and R functions do not match")
all(isapprox.(x2[:ev_Lc],ev_Lc_2)) || error("julia and R functions do not match")
isapprox(y,y2) || error("julia and R functions do not match")

##------------------------------------------------------------------------------
## Projection Kernel: 
##------------------------------------------------------------------------------
R"""
z1 = proj_GAMuT_pheno(Y1)
z2 = proj_GAMuT_pheno(Y1)
w = TestGAMuT(z1$Kc, z1$ev_Kc, z2$Kc, z2$ev_Kc)
"""
@rget Y1
@rget Y2
@rget w

kc,ev_Kc = proj_GAMuT_pheno(Y1)
kc_2,ev_Kc_2 = proj_GAMuT_pheno(Y2)
w2 = testGAMuT(kc,ev_Kc,kc_2,ev_Kc_2)

isapprox(w, w2) || error("julia and R functions do not match")