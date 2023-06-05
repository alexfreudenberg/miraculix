! License: MiXBLUP
!
! Authors 
! Jeremie Vandenplas, jeremie.vandenplas@wur.nl
!
! Copyright (C) 2022-2023 Jeremie Vandenplas
!

module modmiraculix_gpu
 use, intrinsic:: iso_c_binding, only: c_char, c_int, c_long, c_double, c_ptr
 implicit none
 private
 public :: c_sparse2gpu
 public :: c_dcsrtrsv_solve_gpu
 public :: c_free_sparse_gpu


 interface
  !Wrapper for sparse_solve_init.
  !void sparse2gpu(double *V, int *I, int *J, long nnz, long m, long max_ncol,
  !                void **GPU_obj, int *status);
  subroutine c_sparse2gpu(V, I, J, nnz, m, max_ncol, GPU_obj, status)&
              bind(C, name='sparse2gpu')
   import c_int, c_double, c_long, c_ptr
   implicit none
   real(c_double), intent(in) :: V(*)
   integer(c_int), intent(in) :: I(*)
   integer(c_int), intent(in) :: J(*)
   integer(c_long), value, intent(in) :: nnz
   integer(c_long), value, intent(in) :: m
   integer(c_long), value, intent(in) :: max_ncol
   type(c_ptr), intent(out) :: GPU_obj
   integer(c_int), intent(inout) :: status
  end subroutine

  !Wrapper for sparse_solve_compute.
  !void dcsrtrsv_solve_gpu(void *GPU_obj, double *B, int ncol, double *X,
  !                        int *status);
  subroutine c_dcsrtrsv_solve_gpu(GPU_obj, B, ncol, X, status)&
              bind(C, name='dcsrtrsv_solve_gpu')
   import c_int, c_double, c_ptr
   implicit none
   type(c_ptr), value, intent(in) :: GPU_obj
   real(c_double), intent(in) :: B(ncol, *)
   integer(c_int), value, intent(in) :: ncol
   real(c_double), intent(in) :: X(ncol, *)
   integer(c_int), intent(inout) :: status
  end subroutine

  !Wrapper for sparse_solve_destroy.
  !void free_sparse_gpu(void **GPU_obj, int *status);
  subroutine c_free_sparse_gpu(GPU_obj, status)&
              bind(C, name='free_sparse_gpu')
   import c_int, c_ptr
   implicit none
   type(c_ptr), intent(in) :: GPU_obj
   integer(c_int), intent(inout) :: status
  end subroutine

  !Wrapper for dense_solve.
  !void potrs_solve_gpu(double *A, unsigned int input_size, double *B,
  !                     unsigned int rhs_cols, double *X, double *logdet,
  !                     int oversubscribe, int *status);


 end interface

end module
