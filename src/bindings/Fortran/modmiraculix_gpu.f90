! License: MiXBLUP
!
! Authors 
! Jeremie Vandenplas, jeremie.vandenplas@wur.nl
!
! Copyright (C) 2022-2023 Jeremie Vandenplas
!

module modmiraculix_gpu
 use, intrinsic:: iso_c_binding, only: c_char, c_int, c_long, c_double, c_ptr&
                                       , c_null_char
 !$ use omp_lib
 implicit none
 private
 public :: c_sparse2gpu
 public :: c_dcsrtrsv_solve_gpu
 public :: c_free_sparse_gpu
 public :: c_solve_gpu


 interface
  !Wrapper for sparse_solve_init.
  !void sparse2gpu(double *V, long *I, long *J, long nnz, long m, long max_ncol,
  !                int is_lower, void **GPU_obj, int *status);
  subroutine c_sparse2gpu(V, I, J, nnz, m, max_ncol, is_lower, GPU_obj, status)&
              bind(C, name='sparse2gpu')
   import c_int, c_double, c_long, c_ptr
   implicit none
   real(c_double), intent(in) :: V(*)
   integer(c_long), intent(in) :: I(*)
   integer(c_long), intent(in) :: J(*)
   integer(c_long), value, intent(in) :: nnz
   integer(c_long), value, intent(in) :: m
   integer(c_long), value, intent(in) :: max_ncol
   integer(c_int), value, intent(in) :: is_lower
   type(c_ptr), intent(out) :: GPU_obj
   integer(c_int), intent(inout) :: status
  end subroutine

  !Wrapper for sparse_solve_compute.
  !void dcsrtrsv_solve_gpu(void *GPU_obj, char transA, double *B, long ncol, double *X,
  !                        int *status);
  subroutine c_dcsrtrsv_solve_gpu(GPU_obj, transA, B, ncol, X, status)&
              bind(C, name='dcsrtrsv_solve_gpu')
   import c_int, c_long, c_char, c_double, c_ptr
   implicit none
   type(c_ptr), value, intent(in) :: GPU_obj
   character(len=1, kind=c_char), value, intent(in) :: transA
   real(c_double), intent(in) :: B(*)
   integer(c_long), value, intent(in) :: ncol
   real(c_double), intent(inout) :: X(*)
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

 interface c_solve_gpu
  module procedure c_solve_gpu_noperm
  module procedure c_solve_gpu_perm
 end interface

contains

subroutine c_solve_gpu_perm(GPU_obj, is_lower, perm, B, ncol, X, status)
 type(c_ptr), intent(in) :: GPU_obj
 logical, intent(in) :: is_lower
 integer(c_int), intent(in) :: perm(:) !Ap(i,:)=A(perm(i),:)
 real(c_double), intent(in) :: B(:,:)
 integer(c_long), intent(in) :: ncol
 real(c_double), intent(out) :: X(:,:)
 integer(c_int), intent(out) :: status

 integer :: i, j
 real(c_double), allocatable :: X_(:,:)

! allocate(X_, source = B(perm, :))

 allocate(X_, mold=B)
 !$omp parallel default(none) shared(X, B, perm) private(i, j)
 !$omp do collapse(2)
 do j = 1, size(B, 2)
  do i = 1, size(B, 1)
   X(i, j) = B(perm(i), j)
  enddo
 enddo
 !$omp end do
 !$omp end parallel

 status = 0

 if(is_lower)then
  !Solve L*L'*X = B
  call c_dcsrtrsv_solve_gpu(GPU_obj, 'n'//c_null_char, X_, ncol, X, status)
  if(status.ne.0)return

  call c_dcsrtrsv_solve_gpu(GPU_obj, 't'//c_null_char, X, ncol, X_, status)
  if(status.ne.0)return
 else
  !Solve U'*U*X = B
  call c_dcsrtrsv_solve_gpu(GPU_obj, 't'//c_null_char, X_, ncol, X, status)
  if(status.ne.0)return

  call c_dcsrtrsv_solve_gpu(GPU_obj, 'n'//c_null_char, X, ncol, X_, status)
  if(status.ne.0)return
 endif

 X(perm, :) = X_

end subroutine

subroutine c_solve_gpu_noperm(GPU_obj, is_lower, B, ncol, X, status)
 type(c_ptr), intent(in) :: GPU_obj
 logical, intent(in) :: is_lower
 real(c_double), intent(in) :: B(:,:)
 integer(c_long), intent(in) :: ncol
 real(c_double), intent(out) :: X(:,:)
 integer(c_int), intent(out) :: status
 
 real(c_double), allocatable :: X_(:,:)

 allocate(X_, mold = X)

 status = 0

 if(is_lower)then
  !Solve L*L'*X = B
  call c_dcsrtrsv_solve_gpu(GPU_obj, 'n'//c_null_char, B, ncol, X_, status)
  if(status.ne.0)return

  call c_dcsrtrsv_solve_gpu(GPU_obj, 't'//c_null_char, X_, ncol, X, status)
  if(status.ne.0)return
 else
  !Solve U'*U*X = B
  call c_dcsrtrsv_solve_gpu(GPU_obj, 't'//c_null_char, B, ncol, X_, status)
  if(status.ne.0)return

  call c_dcsrtrsv_solve_gpu(GPU_obj, 'n'//c_null_char, X_, ncol, X, status)
  if(status.ne.0)return
 endif

end subroutine

end module
