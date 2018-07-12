## Include necessary functions:
include("../src/functions.jl")

##------------------------------------------------------------------------------
## Variables:
##------------------------------------------------------------------------------
maf = 0.25
traitcor = "med"
nassoc1 = 1
nassoc2 = 1
npheno1 = 1
npheno2 = 1
n_unrelated = 5000
n_variants = 500
causal_var = 0.01
test_approach = 1
sim_approach = 3
ignoreZ = true
srand(1234)


##------------------------------------------------------------------------------
## Simulating genotypes:
##------------------------------------------------------------------------------

UNR_OBS = [size(n_unrelated,1),repeat(NaN,outer=n_variants*n_unrelated)]
for i in 1:n_variants
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    UNR_OBS[:i]= [out[G_mom]]
end
n_causal = floor(causal_var*n_variants)
causal_ind = sample(1:n_variants,n_causal,replace=false)

## determine the MAF of each variant in the sample
MAF_unr = mean(UNR_OBS, 2)/2

## we set "variant" to use the same functions as with cosi
if maf > 0.05
    variant = "common"
else
    variant = "rare"
end


##------------------------------------------------------------------------------
## Simulating phenotypes:
##------------------------------------------------------------------------------
out2 = simulatePhenotypesMediation(npheno1, npheno2,traitcor, nassoc1, nassoc2, causal_ind, MAF_unr, n_unrelated, variant, UNR_OBS,approach=sim_approach, ignoreZ=ignoreZ)
Y1 = out2[Y1]
Y2 = out2[Y2]
