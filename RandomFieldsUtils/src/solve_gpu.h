

#ifndef RFutils_gpusolve
#define RFutils_gpusolve 1


int cholGPU(bool copy, double *matrix, Uint size, double *B, Uint rhs_cols,
     double *LogDet, double *RESULT);
void mgpuSolve(double *matrix, Uint individuals, double *vector);
void gpu_relmat_custom(Uint*, double*, Uint, Uint);
void gpu_relmat_cublas(Uint*, double*, Uint, Uint);

#define PADDIM 4L
#define THREADS_PER_BLOCK 1024 //2048 / 32
#define BLOCKS 1024
#endif
