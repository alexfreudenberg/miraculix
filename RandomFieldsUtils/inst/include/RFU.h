#ifndef RFU_rfutils_h
#define RFU_rfutils_h 1

#include "Basic_utils.h"
#include "errors_messages.h"


#define SCALAR_RU_H 1
#define SCALAR_NEARFMA 5  // never change number, see haplogeno.R !!
#define SCALAR_AVX 6
#define SCALAR_KAHAN 8
#define SCALAR_BASE 2
#define SCALAR_AVX_PARALLEL 9
#define SCALAR_BASE_PARALLEL 10


#define SOLVE 0
#define MATRIXSQRT 1
#define DETERMINANT 2

typedef enum InversionMethod { 
  Cholesky, // 0
  SVD,  // 1
  Eigen, // 2
  Sparse, // 3
  NoInversionMethod, // 4, last user available method
  QR, // 5
  LU, // 6 currently not propagated
  NoFurtherInversionMethod, // 7, local values
  GPUcholesky,		    // 8
  Rcholesky,                // 9
  direct_formula,           // 10
  Diagonal // 10, always last one!
} InversionMethod;
#define nr_InversionMethods ((int) Diagonal + 1)
#define nr_user_InversionMethods ((int) NoFurtherInversionMethod + 1)


#define SOLVE_METHODS 3
typedef // benoetigt
struct solve_storage {
  errorstring_type err_msg;
  InversionMethod method, newMethods[SOLVE_METHODS];
  usr_bool sparse;
  int
  //   SICH_n, MM_n, VT_n, U_n, D_n, 
  //    work_n, w2_n, lnz_n, w3_n, result_n,  nsuper, nnzlindx,
    size, actual_size, actual_pivot,
    nsuper,
    main_n, rhs_n, w2_n, U_n, D_n, w3_n, lnz_n, result_n, 
    
    *pivot_idx, pivot_idx_n, 
    *iwork, iwork_n, //eigen, svd, LU, spam
    *pivotsparse, pivotsparse_n, *xlnz,xlnz_n, //spam
    *snode,snode_n, *xsuper, xsuper_n,*invp,invp_n,   // spam
    *cols,cols_n, *rows,rows_n, *lindx, lindx_n, // spam
    *xja,  xja_n; // chol, eigen, spam
  double 
  *main, *rhs,// diagonal, general -- FORBIDDEN for further use
    *w2, // eigen, svd, LU, QR, pivot
    *U, // eigen, svd, pivot
    *D, // eigen, svd, cholesky, spam, pivot
    *w3, // spam, QR, svd, eigen
    *lnz, // spam, svd
    *result,  // sqrtPosDefFree
    *to_be_deleted; 
} solve_storage;



#endif
