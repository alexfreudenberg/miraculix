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
