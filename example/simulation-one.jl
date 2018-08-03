## Julia script to run one simulation of GAMuT test
## It requires the input variables to be specified in "variables.jl"
## Anna Voss, Claudia Solis-Lemus, July 2018

## Modified to run in HGCC.
## We need to have in the same folder:
## - simulation-one.jl
## - functions.jl
## - variables.jl
## - all_gamut_functions.r

## Replicate to run in case we want to run multiple replicates
## with the same setting:
irep = 1
if(!isempty(ARGS))
    irep = ARGS[1]
end

## Include necessary functions:
include("functions.jl")

##------------------------------------------------------------------------------
## Variables:
##------------------------------------------------------------------------------
include("variables.jl")
using Dates
seed = 3*irep * Dates.hour(now())*Dates.minute(now())*Dates.millisecond(now())
srand(seed);

##------------------------------------------------------------------------------
## Simulating genotypes:
##------------------------------------------------------------------------------
G = zeros(n_unrelated,n_variants)
for i in 1:n_variants
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    G[:,i]= out[2]
end
n_causal = floor(causal_var*n_variants)
causal_ind = sample(collect(1:n_variants),Int(n_causal), replace=false)

## determine the MAF of each variant in the sample
MAF = mapslices(mean, G, 2)/2

##for comparison with Mike's code
using StatsFuns
beta_weight = betapdf(MAF, 1, 25)/betapdf(0, 1, 25)

## we set "variant" to use the same functions as with cosi
variant = maf > 0.05 ? "common" : "rare"

##------------------------------------------------------------------------------
## Simulating phenotypes:
##------------------------------------------------------------------------------
Y = simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF, n_unrelated, G, effectSize)


##------------------------------------------------------------------------------
## GAMuT test
##------------------------------------------------------------------------------
lc,ev_Lc = linear_GAMuT_geno(G)
lc_2,ev_Lc_2 = linear_GAMuT_geno(Y)
pval = testGAMuT(lc,ev_Lc,lc_2,ev_Lc_2)


@rput G
@rput Y
R"""
source("all_gamut_functions.r")
x1 = linear_GAMuT_geno(G)
x2 = linear_GAMuT_geno(Y)
pvalR = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""
@rget pvalR


## Output file
outname = string("gamut-",irep,".txt")
df = DataFrame(i=irep, seed=seed, maf=maf, traitcor=traitcor, nassoc=nassoc, npheno=npheno, n=n_unrelated, nvar = n_variants, causalvar = causal_var, pvalJulia = pval, pvalR = pvalR)

colnames = irep == 1 ? true : false

using CSV
CSV.write(outname, df, header=colnames);
