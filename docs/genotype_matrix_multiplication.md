# Genotype Matrix Multiplication 
The genotype matrix multiplication functionality is exposed through the [5codesAPI](../src/miraculix/5codesAPI.c). Here, we provide a brief overview of the interface:

```{C}
// Options function steering the execution of the matrix multiplication
void setOptions_compressed(
    int use_gpu, // 0 == use CPU
    int cores,   // 0 == use global env
    int floatLoop, // 0 == use doubles
    int meanSubstract, // 1 == increase numerical precis.
    int ignore_missings,// 0==respect missings (untested)
    int do_not_center, // 0 == do centering
    int do_normalize,  // 0 == do not normalize
    int use_miraculix_freq, // 0==use Fortran freq 
    int variant,  // 0 == use internal {0,32,128,256}
    int print_details //1 == get some messages
);

// Preprocess the genotype data and store it
void plink2compressed(
    char *plink, // Pointer to the loaded PLINK bed data (without header bytes)
    char *plink_transposed, // Pointer to the transposed PLINK data 
	int snps, // Number of SNPs
    int indiv, // Number of individuals
    double *f, // Pointer to allele frequencies -- can be omitted with the options above
    int max_n, // Maximal number of columns of matrix to be multiplied
	void **compressed // Pointer to object storing the preprocessed data 
);

// Compute the genotype matrix multiplication
void dgemm_compressed(
    char *trans, // 'N' if A * B or 'T' if A^T * B  
    void *compressed, // Object which has been prepared by the routines above
    int n, // Number of columns of the matrix B
    double *B, // Pointer to matrix B of size snps x n if trans='N' or indiv x n if trans='T'	
    int Ldb, // Leading dimension of B
    double *C, // Result matrix  of size indiv x n if trans='N' or snps x n if trans='T'
    int Ldc // Leading dimension of C
);

// Free the allocated memory
void free_compressed( 
    void **compressed
)

```
Exemplary usage of these interfaces can be found in the Fortran [benchmarking script](../utils/benchmark/benchmark.f90). To get an idea of the concepts, a summary of key parameters is outlined below.

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

