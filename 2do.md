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

1) The c_davies.jl file runs an infinite loop when I run it. No errors are produced, and I know exactly where the problem is, but I want to get your opinion on how to go about changing the code so that it is still relatively similar to the C++ code. 
 
2) The simulationPhenotypesMediation.jl in the seperate.so.tests folder may need to be reviewed. I had made a note earlier that there may be a R function which is not present in Julia. I will look into this more this afternoon. 
 
3) I am currently running simulation-one, so let’s hope that I do not get any errors. I ran it on Friday, but Julia functions were not installed correctly so I had to rerun it a few times.
 
