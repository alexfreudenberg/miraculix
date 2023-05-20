program test_sparse_plink
 use, intrinsic :: iso_c_binding, only: c_loc, c_int, c_double, c_ptr, c_null_ptr
 use, intrinsic :: iso_fortran_env, only: int8, real64
 use modtestplink, only: tgeno
 use mod5codesapi, only: c_sparse_times_plink
 implicit none
 integer(c_int), parameter :: a12_ia(*) = [1, 5, 8]
 integer(c_int), parameter :: a12_ja(*) = [1, 2, 3, 5, 1, 2, 5]
 real(c_double), parameter :: a12_a(*) = [0.5_c_double, 0.5_c_double&
                                          , -1.0_c_double, 0._c_double&
                                          , -1.0_c_double, 0.5_c_double&
                                          , -1._c_double]  
 integer(c_int), parameter :: n1 = size(a12_ia)-1
 integer(c_int), parameter :: n2 = maxval(a12_ja)

 logical, parameter :: ltest = .true.
 logical, parameter :: lcenter = .false.

 integer(c_int) :: nsnp_, nan_
 integer(c_int) :: nIdx
 integer(c_int) :: ldc
 integer(kind=int8), pointer, contiguous :: ivalt(:,:)
 real(kind=c_double), allocatable :: C(:,:)
 real(kind=real64), allocatable :: freq(:)
 character(len=:), allocatable :: freqfile
 type(tgeno) :: mat
 type(tgeno) :: mat_ascii
 type(c_ptr) :: c_plink


 freqfile = 'geno.freq'
 mat%genfile = 'geno.bed'

 write(*,'(/2a)')' Bed file        : ',trim(mat%genfile)
 write(*,'(2a/)')' Allele freq file: ',trim(freqfile)

 call mat%read_bed()
 freq = readfreq(freqfile, mat%nsnp)

 write(*,'(a,i0)')'Number of non-genotyped animals (row sparse): ', n1
 write(*,'(a,i0)')'Number of genotyped animals (column sparse) : ', n2

 if(ltest)then
  mat_ascii%genfile = mat%genfile
  call mat_ascii%read_asci()
  if(lcenter)call center_by_f(mat_ascii, freq)

  associate(densesp => sp_to_dense(n1, n2, a12_ia, a12_ja, a12_a))
   write(*, '(a)')'TEST: Sparse matrix'
   call printdense(densesp)
   write(*, '(a)')'TEST: Sparse matrix * genotype'
   call printdense(matmul(densesp, mat_ascii%cov%val))
  end associate
 endif

 !miraculix
 nan_ = mat%cov%nrows
 nsnp_ = mat%cov%neffectivelevel
 write(*,'(a,i0)')'Number of genotyped animals (plink) : ', nan_
 write(*,'(a,i0)')'Number of SNPs (plink)              : ', nsnp_

 nIdx = n1
 ldc = n1
 write(*,'(a,i0)')'nIdx: ', nIdx
 write(*,'(a,i0)')'ldc : ', ldc

 allocate(ivalt, source = transpose(mat%cov%ival))
 c_plink = c_loc(ivalt)

 call c_sparse_times_plink('N', 'N', c_plink, c_null_ptr&
                           , nsnp_, nan_, nIdx&
                           , a12_ia, a12_ja, a12_a&
                           , C, ldc)


contains

function readfreq(freqfile, nsnp) result(freq)
 character(*), intent(in) :: freqfile
 integer, intent(in) :: nsnp
 real(real64) :: freq(nsnp)

 integer :: i, j, un

 open(newunit=un, file=trim(freqfile), status='old', action='read')
 do i = 1, nsnp
  read(un, *)j, freq(i)
 enddo
 close(un)

end function

pure subroutine center_by_f(mat, freq)
 type(tgeno), intent(inout) :: mat
 real(real64), intent(in) :: freq(:)

 integer :: i, j

 do i = 1, mat%cov%ncol
  do j = 1, mat%cov%nrows
   if(mat%cov%val(j,i).le.2._real64+epsilon(1._real64))then
    mat%cov%val(j,i) = mat%cov%val(j,i) - 2._real64 * freq(j)
   endif
  enddo
 enddo

end subroutine

function sp_to_dense(n1, n2, ia, ja, a) result(mat)
 integer(c_int), intent(in) :: n1
 integer(c_int), intent(in) :: n2
 integer(c_int), intent(in) :: ia(:)
 integer(c_int), intent(in) :: ja(:)
 real(c_double), intent(in) :: a(:)
 real(real64) :: mat(n1, n2)

 integer :: i, j

 mat = 0
 do i = 1, size(ia)-1
  do j = ia(i), ia(i+1)-1
   mat(i, ja(j)) = a(j)
  enddo
 enddo

end function

subroutine printdense(mat)
 real(real64), intent(in) :: mat(:,:)

 integer :: i

 write(*, '(a)')repeat('*', 5)
 do i = 1, size(mat, 1)
  write(*, '(*(f0.3, 1x))')mat(i,:)
 enddo
 write(*, '(a)')repeat('*', 5)

end subroutine

end program
