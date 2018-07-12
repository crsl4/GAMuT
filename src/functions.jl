## Dependencies:
using Distributions
using MultivariateStats
using GLM
using NamedArrays
using DataFrames
using RCall
using Rmath
using StatsBase


## Including necessary functions:
include("seperate.so.tests/Parameters4PhenotypeSimulation.jl")
include("seperate.so.tests/simulateFullyMediatedPhenotype.jl")
include("seperate.so.tests/blockMatrixAntiDiagonal.jl")
include("seperate.so.tests/simulatePartiallyMediatedPhenotype.jl")
include("seperate.so.tests/blockMatrixDiagonal.jl")
include("seperate.so.tests/simulatePhenotypes.jl")
include("seperate.so.tests/createCovMatrix.jl")
include("seperate.so.tests/simulatePhenotypesMediation.jl")
include("seperate.so.tests/minorAlleleCountsBinomial.jl")
include("seperate.so.tests/splitPhenotypeMatrix.jl")

include("linear_GAMuT_geno.jl")
include("proj_GAMuT_pheno.jl")
include("test_GAMuT.jl")
