program test_solve
 use, intrinsic:: iso_c_binding, only: c_int, c_long, c_double, c_ptr
 use modmiraculix_gpu, only: c_sparse2gpu, c_dcsrtrsv_solve_gpu&
                             , c_free_sparse_gpu
 implicit none
 integer(c_int), parameter :: I(*) = [1_c_int, 1_c_int, 2_c_int]
 integer(c_int), parameter :: J(*) = [1_c_int, 2_c_int, 2_c_int]
 real(c_double), parameter :: V(*) = [sqrt(2._c_double)&
                                      , -(sqrt(2._c_double)/2._c_double)&
                                      , sqrt(2._c_double)]
 integer(c_long), parameter :: nnz = size(V)
 integer(c_long), parameter :: m = int(max(maxval(I), maxval(J)), c_long)
 integer(c_long), parameter :: max_ncol = 12
 integer(c_int), parameter :: ncol = 12

 integer(c_int) :: status
 real(c_double) :: B(ncol, m) 
 real(c_double) :: X(ncol, m)
 type(c_ptr) :: GPU_obj


 call c_sparse2gpu(V, I, J, nnz, m, max_ncol, GPU_obj, status)

 call random_number(B)

 call c_dcsrtrsv_solve_gpu(GPU_obj, B, ncol, X, status)

 call printdense(B)
 call printdense(X)

 call c_free_sparse_gpu(GPU_obj, status)

contains

subroutine printdense(mat)
 real(c_double), intent(in) :: mat(:,:)

 integer :: i

 write(*, '(a)')repeat('*', 5)
 do i = 1, size(mat, 1)
  write(*, '(*(f0.3, 1x))')mat(i,:)
 enddo
 write(*, '(a)')repeat('*', 5)

end subroutine

end program
