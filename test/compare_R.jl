using RCall
##------------------------------------------------------------------------------
## Linear Kernel: 
##------------------------------------------------------------------------------
R"""
source("./src/all_gamut_functions.r")

a = matrix(c(0.182378, 0.764705, 0.351618, 0.599969, 0.786265, 0.0844763, .043153, 0.969674, 0.408121, 0.233296, 0.250268, 0.641506), nrow=4, ncol=3, byrow=TRUE)
b = matrix(c(0.764994, 0.906526, 0.0901527, 0.00896168, 0.944189, 0.151536, 0.398221, 0.445275, 0.0962149, 0.994083, 0.969015, 0.976275), nrow=4, ncol=3, byrow=TRUE)

x1 = linear_GAMuT_geno(a)
x2 = linear_GAMuT_geno(b)
y = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""
@rget a
@rget b
@rget x1
@rget x2
@rget y 

include("./src/linear_GAMuT_geno.jl")
include("./src/proj_GAMuT_pheno.jl")
include("./src/test_GAMuT.jl")
lc,ev_Lc = linear_GAMuT_geno(a)
lc_2,ev_Lc_2 = linear_GAMuT_geno(b)
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
z1 = proj_GAMuT_pheno(a)
z2 = proj_GAMuT_pheno(b)
w = TestGAMuT(z1$Kc, z1$ev_Kc, z2$Kc, z2$ev_Kc)
"""

@rget w
kc,ev_Kc = proj_GAMuT_pheno(a)
kc_2,ev_Kc_2 = proj_GAMuT_pheno(b)
w2 = testGAMuT(kc,ev_Kc,kc_2,ev_Kc_2)

isapprox(w, w2) || error("julia and R functions do not match")