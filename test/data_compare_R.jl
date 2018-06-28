using RCall
using Base.Test
##------------------------------------------------------------------------------
## Linear Kernel: 
##------------------------------------------------------------------------------
R"""
source("./src/all_gamut_functions.r")
source("./src/simulation-one.r")

##x1 = linear_GAMuT_geno(Y1)
##x2 = linear_GAMuT_geno(Y2)
x1 = linear_GAMuT_geno(X)
x2 = linear_GAMuT_geno(Y)
y = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""
##@rget Y1
##@rget Y2
@rget X
@rget Y
@rget x1
@rget x2
@rget y 

include("./src/linear_GAMuT_geno.jl")
include("./src/proj_GAMuT_pheno.jl")
include("./src/test_GAMuT.jl")

##lc,ev_Lc = linear_GAMuT_geno(Y1)
##lc_2,ev_Lc_2 = linear_GAMuT_geno(Y2)
lc,ev_Lc = linear_GAMuT_geno(X)
lc_2,ev_Lc_2 = linear_GAMuT_geno(Y)
y2 = testGAMuT(lc,ev_Lc,lc_2,ev_Lc_2)



##Test
@testset "Linear Kernel" begin
    @test x1[:Lc] == lc
    @test x1[:ev_Lc] ≈ ev_Lc[1]
    @test x2[:Lc] == lc_2
    @test x2[:ev_Lc] ≈ ev_Lc_2[1]
    @test y == y2
end 

##------------------------------------------------------------------------------
## Projection Kernel: 
##------------------------------------------------------------------------------
R"""
z1 = proj_GAMuT_pheno(Y1)
z2 = proj_GAMuT_pheno(Y2)
w = TestGAMuT(z1$Kc, z1$ev_Kc, z2$Kc, z2$ev_Kc)
"""
##@rget Y1
##@rget Y2
@rget z1
@rget z2
@rget w

kc,ev_Kc = proj_GAMuT_pheno(Y1)
kc_2,ev_Kc_2 = proj_GAMuT_pheno(Y2)
w2 = testGAMuT(kc,ev_Kc,kc_2,ev_Kc_2)

##Test
@testset "Projection Kernel" begin
    @test z1[:Kc] ≈ kc  ## weird error with ==
    @test z1[:ev_Kc] ≈ ev_Kc[1]
    @test z2[:Kc] ≈ kc_2
    @test z2[:ev_Kc] ≈ ev_Kc_2[1]
    @test w == w2
end 
