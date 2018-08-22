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
using Random
seed = 3*parse(Int64, irep)*hour(now())*minute(now())*millisecond(now())
Random.seed!(seed);

##------------------------------------------------------------------------------
## Simulating genotypes:
##------------------------------------------------------------------------------
G = zeros(n_unrelated,n_variants)
for i in 1:n_variants
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    G[:,i]= out[2]
end

#Not compatible with 0.7
## determine the MAF of each variant in the sample
#MAF = mapslices(mean, G, dims=1)/2

##for comparison with Mike's code
#beta_weight = dbeta.(MAF, 1, 25)/dbeta.(0, 1, 25)
#G0 = G * diagm(vec(beta_weight))
#G = G0 .- mean(G0, 1)

##------------------------------------------------------------------------------
## Simulating phenotypes:
##------------------------------------------------------------------------------
n_causal = floor(causal_var*n_variants)
causal_ind = sample(collect(1:n_variants),Int(n_causal), replace=false)

## we set "variant" to use the same functions as with cosi
variant = maf > 0.05 ? "common" : "rare"

Y = simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF, n_unrelated, G, effectSize)

# Not compatible with 0.7
#P0 = scale!(Y, 2)
#P = P0 .-mean(P0,1)
##------------------------------------------------------------------------------
## GAMuT test
##------------------------------------------------------------------------------
lc,ev_Lc = linear_GAMuT_geno(Y)
lc_2,ev_Lc_2 = linear_GAMuT_geno(G) ##produces only 2 values rather than a large matrix when it works 
pval = testGAMuT(lc,ev_Lc,lc_2,ev_Lc_2)


@rput G
@rput Y
R"""
source("./src/r-scripts/all_gamut_functions.r")
x1 = linear_GAMuT_geno(Y)
x2 = linear_GAMuT_geno(G)
pvalR = TestGAMuT(x1$Lc,x1$ev_Lc,x2$Lc,x2$ev_Lc)
"""
@rget pvalR


## Output file
outname = string("gamut-",irep,".txt")
df = DataFrame(i=irep, seed=seed, maf=maf, traitcor=traitcor, nassoc=nassoc, npheno=npheno, n=n_unrelated, nvar = n_variants, causalvar = causal_var, pvalJulia = pval, pvalR = pvalR)

colnames = irep == 1 ? true : false

#Not compatible with 0.7
#using CSV
#CSV.write(outname, df, header=colnames);

g = open("file.txt","w")
for i in 1:nrows
    write(g,string(df,"\n"))
end
close(g)