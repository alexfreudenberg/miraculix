/*
 Authors 
 Alexander Freudenberg alexander.freudenberg@uni-mannheim.de

 Copyright (C) 2022-2023 Alexander Freudenberg

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

void plink2gpu(
    char *plink,            // Pointer to SNP data in PLINK format
    char *plink_transposed, // Pointer to transposed SNP data in PLINK format
    int snps,               // Number of SNPs
    int indiv,              // Number of individuals
    double *f,              // Vector of allele frequencies
    int max_ncols,          // Maximum number of columns of RHS
    void **GPU_obj);        // Pointer to struct with GPU data info
void dgemm_compressed_gpu(
    bool trans,     // If true, transposed matrix multiplication will be used
    void *GPU_obj,  // Pointer to struct with GPU data info
    int n,          // Actual number of columns of matrix B
    double *B,      // RHS matrix of dimensions (nsnps, n) if trans =
                    // false and (nindiv,n) if trans = true
    int ldb,        // Not used
    int centered,   // SNP data will be centered if centerd = 1
    int normalized, // Not used currently
    double *C,      // Output matrix
    int ldc         // Not used
);
void freegpu(void **GPU_obj);

