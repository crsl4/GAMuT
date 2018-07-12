## Test script to compare julia code and R code for GAMuT test
## It is not an automatic test, as we need to specify to travis
## that R needs to have CompQuadForm package installed (and we
## do not know how to do this yet).
## Anna Voss, Claudia Solis-Lemus, July 2018

using Base.Test

## Include necessary functions:
include("../src/functions.jl")

## Input variables:
maf = 0.25
traitcor = "med"
nassoc = 1
npheno = 1
n_unrelated = 5000
n_variants = 500
causal_var = 0.01

@rput maf
@rput traitcor
@rput nassoc
@rput npheno
@rput n_unrelated
@rput n_variants
@rput causal_var


##------------------------------------------------------------------------------
## First, we create phenotype and genotype matrices in R, and copy back to Julia;
##------------------------------------------------------------------------------

## Linear Kernel:
R"""
set.seed(1234)
source("../src/r-scripts/all_gamut_functions.r")
source("../src/r-scripts/simulation-one.r")

x1 = linear_GAMuT_geno(X)
x2 = linear_GAMuT_geno(Y)
pvalR = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""

@rget X
@rget Y
@rget x1
@rget x2
@rget pvalR


## Now running in Julia:
lc,ev_Lc = linear_GAMuT_geno(X)
lc_2,ev_Lc_2 = linear_GAMuT_geno(Y)
pvalJulia = testGAMuT(lc,ev_Lc,lc_2,ev_Lc_2)

##Test
@testset "Linear Kernel from R to Julia" begin
    @test x1[:Lc] ≈  lc
    @test x1[:ev_Lc] ≈ ev_Lc
    @test x2[:Lc] ≈ lc_2
    @test x2[:ev_Lc] ≈ ev_Lc_2[1]
    @test pvalR ≈ pvalJulia
end

##------------------------------------------------------------------------------
## Second, we will create all matrices in Julia and copy to R to compare performance
##------------------------------------------------------------------------------
srand(1234)

X = zeros(n_unrelated,n_variants)
for i in 1:n_variants
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    X[:,i]= out[2]
end
n_causal = floor(causal_var*n_variants)
causal_ind = sample(collect(1:n_variants),Int(n_causal), replace=false)

## determine the MAF of each variant in the sample
MAF = mean(X, 2)/2

## we set "variant" to use the same functions as with cosi
variant = maf > 0.05 ? "common" : "rare"

## Simulating phenotypes:
Y = simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF, n_unrelated, X)

## Test GAMuT in julia:
lc,ev_Lc = linear_GAMuT_geno(X)
lc_2,ev_Lc_2 = linear_GAMuT_geno(Y)
pvalJulia = testGAMuT(lc,ev_Lc,lc_2,ev_Lc_2)

@rput X
@rput Y
R"""
x1 = linear_GAMuT_geno(X)
x2 = linear_GAMuT_geno(Y)
pvalR = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""

@rget x1
@rget x2
@rget pvalR

@testset "Linear Kernel from Julia to R" begin
    @test x1[:Lc] ≈  lc
    @test x1[:ev_Lc] ≈ ev_Lc
    @test x2[:Lc] ≈ lc_2
    @test x2[:ev_Lc] ≈ ev_Lc_2[1]
    @test pvalR ≈ pvalJulia
end

##################################################
## Below not tested yet:

##------------------------------------------------------------------------------
## Projection Kernel:
##------------------------------------------------------------------------------
R"""
z1 = proj_GAMuT_pheno(X)
z2 = proj_GAMuT_pheno(Y)
pvalR2 = TestGAMuT(z1$Kc, z1$ev_Kc, z2$Kc, z2$ev_Kc)
"""

@rget z1
@rget z2
@rget pvalR2

kc,ev_Kc = proj_GAMuT_pheno(X)
kc_2,ev_Kc_2 = proj_GAMuT_pheno(Y)
pvalJulia2 = testGAMuT(kc,ev_Kc,kc_2,ev_Kc_2)

##Test
@testset "Projection Kernel" begin
    @test z1[:Kc] ≈ kc  ## weird error with ==
    @test z1[:ev_Kc] ≈ ev_Kc[1]
    @test z2[:Kc] ≈ kc_2
    @test z2[:ev_Kc] ≈ ev_Kc_2[1]
    @test pvalR2 ≈ pvalJulia2
end
