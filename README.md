# miraculix
#### Efficient algebraic functions for statistical analysis of genomic data
*May 2023*

## Description
miraculix is a C/CUDA library for mathematical operations on compressed genotype data. It provides highly efficient routines that can be used in statistical analyses of genomic information, e.g. genome-wide association studies, genomic breeding value estimation and population summary statistics. As such, it offers interfaces to allow integration into existing C utilities as well as interoperability with higher-level programming languages such as R, Julia or Fortran.

Through its CUDA implementation, miraculix aims to open the high-performance computing capabilities of modern Nvidia GPUs to researchers and practitioners in the field of statistical genomics. To this end, miraculix uses and extends [Nvidia's CUTLASS (2.10)](https://github.com/NVIDIA/cutlass) library to enable efficient data movement on the GPU. [Experiments]() suggest that, for instance, genomic breeding value estimation benefits greatly from offloading bottlenecks to GPUs.

Previous versions of miraculix have been released as an R package on CRAN. However, CRAN's strict portability requirements have been shown to be too restrictive to maintain the code base while simultaneously assuring efficiency across instruction set architectures. Yet, many of the interfaces required for a smooth R integration are still present and can be compiled to a full R package.

While we migrate the existing functionality to this repository, please expect further changes to its structure.

## Requirements
For compilation, the following software is required:
* [GNU Make](https://www.gnu.org/software/make/) *(other implementations of Make have not been tested)*
* GNU's gcc or Intel's icc compiler

Additionally, to compile the CUDA code, you need
* a CUDA installation (11.0 or newer)
  

## Compilation
To compile the CUDA version of miraculix, you first need to clone the CUTLASS submodule:
```{bash}
git submodule init
git submodule update
```
Then, the miraculix library is compiled through
```{bash}
cd src
make libmiraculix
```
The miraculix CUDA implementation is compiled separately through
```{bash}
make libmiraculixGPU
```

## Available Interfaces
If you wish to integrate the genotype matrix multiplication functionality into your genomic analysis pipeline, please use the dedicated interfaces. You can find an overview of them [here](./docs/genotype_matrix_multiplication.md). 

Interfaces to the other routines will be added gradually.

## Examples 
Exemplary usage of the routines can be found in the [examples](./examples) folder and in the [benchmarking files](./utils/benchmark/). The syntax for these Fortran files is as follows:
**Options**
Set the options used later on. Most of the options can't be changed after they have been set initially. 
```{Fortran}
call c_setOptions_compressed(usegpu, 0, 0, 0, 1, c_not_center, 0, 0, c_variant, c_lverbose)
```
Through the `usegpu` parameter, we can control if the GPU implementation or the *5codes* algorithm on the CPU is used. 
The `c_not_center` parameter turns centering of the genotype matrix off, if it is set to 1. 
The `c_variant` parameter chooses which internal implementation is used - this can be experimented with to find the optimal performance. 
Through `c_lverbose` we control how many internal information is printed to stdout.
**Initialization**
Preprocess SNP data and store it in a separate object.
```{Fortran}
call c_plink2compressed(c_plinkbed, c_plinkbed_transposed, c_snps, c_indiv, c_f, c_ncol, c_compressed)
```
The `c_plinkbed` and `c_plinkbed_transposed` hold the SNP data in [PLINK bed](https://www.cog-genomics.org/plink/1.9/formats#bed) format truncated by the header bytes. The number of SNPs and individuals in the dataset is supplied through `c_snps` and `c_indiv`. Allele frequencies can be supplied through the `c_f` pointer - this can be useful when using frequencies that differ from the SNP data. The `c_ncol` parameter is used for the GPU implementation: It indicates the maximum number of columns with which the `dgemm_compressed` routine will be called. The `c_compressed` parameter holds a pointer to a pointer, in which the preprocessed data storage object will be stored. 

**Computation**
Calculate the genotype matrix multiplication.
```{Fortran}
call c_dgemm_compressed('n', c_compressed, c_ncol, c_B, c_snps, c_C, c_indiv)
```
In this example, we calculate the untransposed (hence 'n') genotype matrix times a real-valued matrix in double precision stored in `c_B`.  

**Destroy**
Free allocated memory.
```{Fortran}
call c_free_compressed(c_compressed)
```

## Citation
If you decide to use this repository for your scientific work, please consider citing it.

## Publications
- **The Modular Breeding Program Simulator (MoBPS) allows efficient simulation of complex breeding programs.** 
Pook, T., Reimer, C., Freudenberg, A., Büttgen, L., Geibel, J., Ganesan, A., Ha, T., Schlather, M., Mikkelsen, L., Simianer, H.,
Animal Production Science, 61, 1982-1989, 2021. 
