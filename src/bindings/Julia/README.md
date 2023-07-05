## Julia interfaces to miraculix
This directory includes the skeleton of a Julia package for processing genotype data using functionality in `miraculix`. Additionally, a number of auxiliary utilities are defined. 

To use this package, you need to consider the notes on compatibility in [Project.toml](./Project.toml) and install the dependencies therein. Unfortunately, installation from within Julia through `instantiate` **is not possible.**

The package is structured as follows:

* [`miraculix.jl`](./miraculix.jl): Defines the core of the package and general purpose utilties used in other modules.
* [`dgemm_compressed.jl`](./dgemm_compressed.jl): Defines an interface to the GPU-based multiplication of a SNP matrix with a floating-point matrix.
* [`crossproduct.jl`](./crossproduct.jl): Defines an interface to the GPU-based SNP matrix multiplication.
* [`solve.jl`](./solve.jl): Provides an interface to GPU-based functionality for solving sparse or dense equation systems.
* [`compressed_operations.jl`](./compressed_operations.jl): Introduces a number of useful operations on compressed genotype data.
* [`read_plink.jl`](./read_plink.jl): Provides functionality to read genotype data stored in PLINK binary format and operations working on this data format.