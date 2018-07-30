We want to compare the performance of R and julia.
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
Pkg.add("MultivariateStats") ## we did not have this
include("simulation-one.jl")

signal (11): Segmentation fault
while loading /mnt/icebreaker/data2/home/csolislemus/gamut/functions.jl, in expression starting on line 7
Rf_install at /usr/lib64/R/lib/libR.so (unknown line)
R_initMethodDispatch at /home/vpatel/Software/Test/R-3.4.3/src/library/methods/src/methods_list_dispatch.c:95
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6962
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:747
do_set at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:2585
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6493
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
forcePromise at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:520
FORCE_PROMISE at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4742 [inlined]
getvar at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:4784 [inlined]
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6233
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
bcEval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:6444
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:624
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:747
Rf_evalList at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:2681
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:719
do_if at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1909
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:700
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:700
do_begin at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:2192
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:700
R_execClosure at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:1614
Rf_eval at /home/vpatel/Software/Test/R-3.4.3/src/main/eval.c:747
setup_Rmainloop at /home/vpatel/Software/Test/R-3.4.3/src/main/main.c:931
Rf_initEmbeddedR at /home/vpatel/Software/Test/R-3.4.3/src/unix/Rembedded.c:63
initEmbeddedR at /home/csolislemus/.julia/v0.6/RCall/src/setup.jl:184
__init__ at /home/csolislemus/.julia/v0.6/RCall/src/setup.jl:222
unknown function (ip: 0x2b61ae6189ff)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
jl_apply at /home/centos/buildbot/slave/package_tarball64/build/src/julia.h:1424 [inlined]
jl_module_run_initializer at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:87
jl_init_restored_modules at /home/centos/buildbot/slave/package_tarball64/build/src/dump.c:2443 [inlined]
_jl_restore_incremental at /home/centos/buildbot/slave/package_tarball64/build/src/dump.c:3272
jl_restore_incremental at /home/centos/buildbot/slave/package_tarball64/build/src/dump.c:3292
_include_from_serialized at ./loading.jl:157
_require_from_serialized at ./loading.jl:200
_require at ./loading.jl:457
require at ./loading.jl:398
unknown function (ip: 0x2b61ae607aa2)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
jl_apply at /home/centos/buildbot/slave/package_tarball64/build/src/julia.h:1424 [inlined]
eval_import_path_ at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:403
eval_import_path at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:430 [inlined]
jl_toplevel_eval_flex at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:495
jl_parse_eval_all at /home/centos/buildbot/slave/package_tarball64/build/src/ast.c:873
jl_load at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:616
include_from_node1 at ./loading.jl:569
unknown function (ip: 0x2b618a515b1b)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
include at ./sysimg.jl:14
unknown function (ip: 0x2b618a3ba86b)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
do_call at /home/centos/buildbot/slave/package_tarball64/build/src/interpreter.c:75
eval at /home/centos/buildbot/slave/package_tarball64/build/src/interpreter.c:242
jl_interpret_toplevel_expr at /home/centos/buildbot/slave/package_tarball64/build/src/interpreter.c:34
jl_toplevel_eval_flex at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:577
jl_parse_eval_all at /home/centos/buildbot/slave/package_tarball64/build/src/ast.c:873
jl_load at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:616
include_from_node1 at ./loading.jl:569
unknown function (ip: 0x2b618a515b1b)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
include at ./sysimg.jl:14
unknown function (ip: 0x2b618a3ba86b)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
do_call at /home/centos/buildbot/slave/package_tarball64/build/src/interpreter.c:75
eval at /home/centos/buildbot/slave/package_tarball64/build/src/interpreter.c:242
jl_interpret_toplevel_expr at /home/centos/buildbot/slave/package_tarball64/build/src/interpreter.c:34
jl_toplevel_eval_flex at /home/centos/buildbot/slave/package_tarball64/build/src/toplevel.c:577
jl_toplevel_eval_in at /home/centos/buildbot/slave/package_tarball64/build/src/builtins.c:496
eval at ./boot.jl:235
unknown function (ip: 0x2b618a4df39f)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
eval_user_input at ./REPL.jl:66
unknown function (ip: 0x2b618a54d1cf)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
macro expansion at ./REPL.jl:97 [inlined]
#1 at ./event.jl:73
unknown function (ip: 0x2b61ae546bbf)
jl_call_fptr_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /home/centos/buildbot/slave/package_tarball64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /home/centos/buildbot/slave/package_tarball64/build/src/gf.c:1933
jl_apply at /home/centos/buildbot/slave/package_tarball64/build/src/julia.h:1424 [inlined]
start_task at /home/centos/buildbot/slave/package_tarball64/build/src/task.c:267
unknown function (ip: 0xffffffffffffffff)
Allocations: 47418575 (Pool: 47407220; Big: 11355); GC: 108
Segmentation fault (core dumped)
```

Problem seems to be `RCall`.




Anna: When I tries to add packages to Julia, specifically Distributions, I get the following error: 

```julia
Pkg.add("Distributions")

ERROR: GitError(Code:ENOTFOUND, Class:Repository, Could not find repository from 'METADATA')
Stacktrace:
 [1] macro expansion at ./libgit2/error.jl:99 [inlined]
 [2] Base.LibGit2.GitRepo(::String) at ./libgit2/repository.jl:10
 [3] macro expansion at ./pkg/entry.jl:57 [inlined]
 [4] macro expansion at ./task.jl:302 [inlined]
 [5] add(::String, ::Base.Pkg.Types.VersionSet) at ./pkg/entry.jl:51
 [6] (::Base.Pkg.Dir.##4#7{Array{Any,1},Base.Pkg.Entry.#add,Tuple{String}})() at ./pkg/dir.jl:36
 [7] cd(::Base.Pkg.Dir.##4#7{Array{Any,1},Base.Pkg.Entry.#add,Tuple{String}}, ::String) at ./file.jl:70
 [8] #cd#1(::Array{Any,1}, ::Function, ::Function, ::String, ::Vararg{String,N} where N) at ./pkg/dir.jl:36
 [9] add(::String) at ./pkg/pkg.jl:117
```

Problem with Github connecting to HGCC. 
I tried creating an ssh key pair through GitHub, but the key pair will not load with HGCC. 