## Dependencies:
using Distributions
using MultivariateStats
using GLM
using DataFrames
using RCall
using Rmath
using StatsBase


## Including necessary functions:
include("Parameters4PhenotypeSimulation.jl")
include("simulateFullyMediatedPhenotype.jl")
##include("seperate.so.tests/blockMatrixAntiDiagonal.jl")
##include("seperate.so.tests/blockMatrixDiagonal.jl")
include("simulatePartiallyMediatedPhenotype.jl")
include("simulatePhenotypes.jl")
include("createCovMatrix.jl")
include("simulatePhenotypesMediation.jl")
include("minorAlleleCountsBinomial.jl")
include("splitPhenotypeMatrix.jl")

include("linear_GAMuT_geno.jl")
include("proj_GAMuT_pheno.jl")
include("test_GAMuT.jl")
