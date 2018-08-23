- test data_compare_R as it is: test for different seeds and variable inputs
- example: simulation-one.jl, do ~1000 simulations per setting. Two settings: nassoc=1 and nassoc=0 (for npheno=2)
- do data_compare_R a true test function, add travis stuff:
```
no we don’t need ape for the tests. You can check our travis file, line 28:
https://github.com/crsl4/PhyloNetworks.jl/blob/master/.travis.yml#L28
where there’s a line to ask travis to install ggplot2. This line is commented out (because it took 6 minutes to install ggplot2: too long, so avoided it in the end). But you can use it as a template to make travis install ape.
```
- put simulation of G and Y into two functions (instead of copying and pasting all the lines everytime)
- add docstrings to all functions
- create test functions independent of R: davies function in julia needed

At the end of Anna's project:
 
 1) In davies.jl, the value for "th" is unclear. th is a vector because it is indexed. th is a pointer in the c code, so I made it an input in the julia code for the order and cfe functions. The value of th is determined at the beginning of the qfc function, but I don't know a Julia function for the c code since it deals with allocation of memory. (see comment in c_davies.jl qfc function)

 1) The c_davies.jl file runs an infinite loop when I run it. Comments are added to the top of c_davies.jl. Running l1 with @show still produces an infinite loop. This indicates that the problem occurs before the if statement in l1. I noticed that there is never a call to "goto" l1 so this might be the source of our problem. 

 3) Run Simulation-one in HGCC. (You may have to go into GitHub to get the Julia 0.6 version of simulation-one.jl)

 4) BTW: I updated the REQUIRE file to match 0.7 (I changed Julia version and deleted MultivariateStats and GLM packages)
