# miraculix
This repository contains open source code for the master thesis 'Statistical Methods for Ultra-High Dimensional Data in Genetics' submitted on 31st of March 2021. 

## Requirements
* R 4.0 or higher
* gcc 9.0 or higher
* CUDA 11.0 or higher
* NVIDIA GPU of compute capability 7.5 or higher (Turing architecture or newer)
* A local copy of cutlass (https://github.com/NVIDIA/cutlass)

## Installation
* R CMD INSTALL RandomFieldsUtils --configure-args="USE_GPU=yes USE_AVX=yes"
* R CMD INSTALL miraculix --configure-args="USE_GPU=yes USE_AVX=yes"
* R CMD INSTALL rpud --configure-args="USE_GPU=yes"

## Structure
* miraculix contains all functions for the calculation of the genomic relationship matrix. GPU code can be found in files starting with the 'mma' prefix.
* RandomFieldsUtils contributes a GPU accelerated Cholesky decomposition to miraculix. The function can be found in solve_gpu.cu
* rpud contains proof-of-concepts for the GPU accelerated REML optimization and phenotypic variance decomposition, GPU pointer management via R, handling of large equation systems.
* R files in the main folder provide benchmarks, use-cases and syntax hints. 

## Known bugs:
* `libcusolver.so.10: cannot open shared object file: No such file or directory`: Execute `export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH` or add it to your bash_source`


![MMAGPU benchmarks](https://github.com/alexfreudenberg/master_thesis_software/raw/master/MMAGPU.png)
![Cholesky benchmarks](https://github.com/alexfreudenberg/master_thesis_software/raw/master/Chol.png)
