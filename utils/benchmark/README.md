In this directory, we define a number of utilities for benchmarking functionality in the `miraculix` library. Unfortunately, there is no comprehensive suite at the moment and different calculations are evaluated by separate files. 

### **GPU-based SNP matrix crossproduct**
This functionaliy is evaluated by the following files:

* `benchmark_suite.jl`: This defines a number of benchmarks for comparing `miraculix` with PLINK and a `cuBLAS` impelmentation and wraps them into a suite.
* `run_suite.jl`: Evaluates specific sets within the suite. It takes the execution mode ("PLINK", "miraculix", "cuBLAS") as a command-line argument.
* `cublas_uint8.cu`: Defines a SNP matrix crossproduct based on 8-bit unsigned integers.


### **GPU-based SNP matrix floating-point multiplication**
This functionaliy is evaluated by the following files:

* `benchmark.f90`: This file is compiled to a binary and takes a path to a genotype dataset in PLINK binary format as input argument.


### **General purpose utilities**

* `Makefile`: compiles the CUDA and Fortran files.
* `benchamrk_*.sh`: These are auxiliary slurm scripts with recommended settings which can be used for evaluating the benchmark scripts on an HPC cluster.