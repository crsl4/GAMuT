We want to compare the performace of R and julia.
We will do this with script `simulation-one.jl` which simulates data, and runs gamut in R and in julia. It creates a data frame at the end with the pvalues.

First, we copy the scripts to HGCC:
```shell
cd Documents/github/GAMuT/example
scp simulation-one.jl csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/gamut/
scp ../src/functions.jl csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/gamut/
scp variables.jl csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/gamut/
scp ../src/r-scripts/all_gamut_functions.r csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/gamut/
```
Ideally, we would want to clone the repo in HGCC, but we do not want to accidentally push output files into the repo.

Now, we login to HGCC:
```shell
ssh csolislemus@hgcc.genetics.emory.edu
```

We want to test the script before running as array job:
```shell
qlogin -q i.q
cd gamut
module load R
module load julia
```
Now, insude julia:
```julia
include("simulation-one.jl")
```