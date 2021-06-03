

#ifndef miraculix_gpusolve
#define miraculix_gpusolve 1

#include "MX.h"

// void gpuSolve(double *matrix, Uint individuals, double *result);
// void gpu_relmat_custom(Uint*, double*, Uint, Uint);
void gpu_relmat_cublas(Uint*, double*, Uint, Uint);

#define PADDIM 4L
#define THREADS_PER_BLOCK 1024 //2048 / 32

#endif
