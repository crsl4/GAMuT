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