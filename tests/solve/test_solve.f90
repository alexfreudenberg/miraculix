program test_solve
 use, intrinsic:: iso_c_binding, only: c_int, c_long, c_double, c_ptr&
                                       , c_null_char
 use modmiraculix_gpu, only: c_sparse2gpu, c_dcsrtrsv_solve_gpu&
                             , c_free_sparse_gpu, c_solve_gpu
 implicit none
 integer(c_long), parameter :: I(*) = [1_c_long, 1_c_long, 2_c_long]
 integer(c_long), parameter :: J(*) = [1_c_long, 2_c_long, 2_c_long]
 real(c_double), parameter :: V(*) = [sqrt(2._c_double)&
                                      , -(sqrt(2._c_double)/2._c_double)&
                                      , sqrt(2._c_double)]
 integer(c_int), parameter :: is_lower = 0_c_int
 integer(c_long), parameter :: nnz = size(V)
 integer(c_long), parameter :: m = int(max(maxval(I), maxval(J)), c_long)
 integer(c_long), parameter :: max_ncol = 12
 integer(c_long), parameter :: ncol = 12
 real(c_double), parameter :: thr = 1000 * epsilon(1._c_double)

 logical, parameter :: lprint = .false.

 integer(c_int) :: i_, status
 real(c_double), allocatable :: mat(:,:)
 real(c_double) :: B(m, ncol) 
 real(c_double) :: X(m, ncol)
 real(c_double) :: X_(m, ncol)
 type(c_ptr) :: GPU_obj

 !Sparse to Dense
 mat = to_dense(I, J, V)

 if(lprint)call printdense(mat, 'LHS:')

 !Init GPU
 call c_sparse2gpu(V, I, J, nnz, m, max_ncol, is_lower, GPU_obj, status)

 !Init RHS
 call random_number(B)
 !B=1
 if(lprint)call printdense(B, 'RHS:')


 !Triangle solve
 call c_dcsrtrsv_solve_gpu(GPU_obj, 'n'//c_null_char, B, ncol, X, status)

 if(lprint)call printdense(X, 'U\B')
 call check(matmul(mat, X) - B, 'U\B')


 !Transpose triangle solve
 call c_dcsrtrsv_solve_gpu(GPU_obj, 't'//c_null_char, B, ncol, X, status)

 if(lprint)call printdense(X, 'U^T\B')
 call check(matmul(transpose(mat), X) - B, 'U^T\B')


 !Solve U'*U*X = B
 call c_dcsrtrsv_solve_gpu(GPU_obj, 't'//c_null_char, B, ncol, X_, status)
 call c_dcsrtrsv_solve_gpu(GPU_obj, 'n'//c_null_char, X_, ncol, X, status)

 if(lprint)call printdense(X, 'U^T\U\B')
 call check(matmul(matmul(transpose(mat), mat), X) - B,  'U^T\U\B')

 !Solve U'*U*X
 X = 0
 call c_solve_gpu(GPU_obj, .false., B, ncol, X, status)

 if(lprint)call printdense(X, 'U^T\U\B')
 call check(matmul(matmul(transpose(mat), mat), X) - B,  'U^T\U\B')


 !Free GPU memory
 call c_free_sparse_gpu(GPU_obj, status)


 !!!!!!!!!!!!!!!!!!!!!
 !Sparse to Dense
 mat = to_dense(J, I, V)
 if(lprint)call printdense(mat, 'LHS:')
 

 !Init GPU
 call c_sparse2gpu(V, J, I, nnz, m, max_ncol, 1_c_int, GPU_obj, status)


 !Solve U'*U*X
 call c_solve_gpu(GPU_obj, .true., B, ncol, X, status)

 if(lprint)call printdense(X, 'L^T\L\B')
 call check(matmul(matmul(mat, transpose(mat)), X) - B,  'L^T\L\B')


 !Solve U'*U*X
 call c_solve_gpu(GPU_obj, .true., [(i_, i_=1, size(B,1))], B, ncol, X, status)

 if(lprint)call printdense(X, 'L^T\L\B')
 call check(matmul(matmul(mat, transpose(mat)), X) - B,  'L^T\L\B')


 !Free GPU memory
 call c_free_sparse_gpu(GPU_obj, status)



contains

subroutine printdense(mat, msg)
 real(c_double), intent(in) :: mat(:,:)
 character(*), intent(in), optional :: msg

 integer :: i

 write(*, '(a)')repeat('*', 5)
 if(present(msg))write(*, '(a)')msg
 do i = 1, size(mat, 1)
  write(*, '(*(f0.5, 1x))')mat(i,:)
 enddo
 write(*, '(a)')repeat('*', 5)

end subroutine

subroutine check(mat, msg)
 real(c_double), intent(in) :: mat(:,:)
 character(*), intent(in) :: msg

 if(all(abs(mat) < thr))then
  write(*, '(a, " OK")')msg
 else
  call printdense(mat, msg//" wrong")
 endif

end subroutine

pure function to_dense(I, J, V) result(mat)
 integer(c_long), intent(in) :: I(:), J(:)
 real(c_double), intent(in) :: V(:)
 real(c_double) :: mat(maxval(I), maxval(J))

 integer :: i_, j_

 mat = 0
 do i_ = 1, size(I)
  mat(I(i_), J(i_)) = V(i_)
 enddo

end function

end program
