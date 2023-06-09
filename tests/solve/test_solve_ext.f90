program test_solve_ext
 use, intrinsic:: iso_c_binding, only: c_char, c_ptr, c_int, c_long&
                                      , c_double
 use modmiraculix_gpu, only: c_sparse2gpu, c_solve_gpu&
                             , c_free_sparse_gpu
 implicit none

 integer(c_int), parameter :: is_lower = 0_c_int
 integer(c_long), parameter :: max_ncol = 8
 integer(c_long), parameter :: ncol = 8

 integer :: un, io
 integer(c_int) :: status
 integer(c_long) :: i_, m
 integer(c_int), allocatable :: perm(:)
 integer(c_long), allocatable :: I(:), J(:)
 real(c_double), allocatable :: V(:)
 real(c_double), allocatable :: B(:,:), X(:,:)
 integer(c_long) :: nnz
 type(c_ptr) :: GPU_obj

 character(30) :: cholfile

 call get_command_argument(1, cholfile)
 write(*,'(a)')'Name of file: ', cholfile

 open(newunit=un, file=cholfile, status='old', action='read')
 nnz = 0
 do
  read(un, *, iostat=io)
  if(io.ne.0)exit
  nnz = nnz + 1 
 enddo
 write(*,'(a,i0)')'Number of non-zero    : ',nnz

 rewind(un)
 allocate(I(nnz))
 allocate(J(nnz))
 allocate(V(nnz))

 do i_ = 1, nnz
  read(un, *, iostat=io)I(i_), J(i_), V(i_)
  if(io.ne.0)exit
 enddo
 close(un)


 !Init GPU
 m = max(maxval(I), maxval(J))
 write(*,'(a,i0)')'Number of rows/columns: ',m
 call c_sparse2gpu(V, I, J, nnz, m, max_ncol, is_lower, GPU_obj, status)

 !Solve U'*U*X
 allocate(B(m, ncol))
 allocate(X(m, ncol))

 call random_number(B)

 call c_solve_gpu(GPU_obj, .false., B, ncol, X, status)

 !Free GPU memory
 call c_free_sparse_gpu(GPU_obj, status)

end program
