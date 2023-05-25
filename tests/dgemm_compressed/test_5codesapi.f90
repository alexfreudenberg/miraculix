!  Authors 
!  Jeremie Vandenplas, jeremie.vandenplas@wur.nl
!  Alexander Freudenberg, alexander.freudenberg@stads.de
!  Copyright (C) 2022-2023 Jeremie Vandenplas, Alexander Freudenberg

!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at

!     http://www.apache.org/licenses/LICENSE-2.0

!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.

! A CUDA installation and a GPU of compute capability greater than 7.0 is required


program test_5codesAPI
 use, intrinsic:: iso_fortran_env, only: int8, real64
 use, intrinsic:: iso_c_binding, only: c_loc, c_int, c_ptr, c_null_ptr, c_double&
                                       , c_char, c_null_char
 use mod5codesapi, only: c_plink2compressed, c_dgemm_compressed, c_free_compressed&
                         , c_get_freq => c_get_compressed_freq&
                         , c_setOptions_compressed
 use modtestplink, only: tgeno
 use modplink_miraculix, only: transpose_integermatrix
 !$ use omp_lib
 implicit none

 integer, parameter :: ncol =  10
 integer, parameter :: nrepet = 3
 integer, parameter :: nmin = merge(1, 0, nrepet > 1)
 logical, parameter :: ltest = .true.
 logical, parameter :: lcenter = .true.
 logical, parameter :: lverbose = .true.


 integer(c_int), parameter :: usegpu = 1_c_int
 integer(c_int), parameter :: c_not_center = merge(0, 1, lcenter)
 integer(c_int), parameter :: c_lverbose = 2!merge(1, 0, lverbose)
 real(real64), parameter :: tol_real64 = 1e-4 !1000 * epsilon(1._real64)
 
 integer :: irepet
 integer :: i, j, io, err, nrow

 integer :: idummy, compressedn, lda_atrans
 integer(kind=int8), allocatable :: trans_ival(:,:)
 character(len=30) :: freqfile
 real(kind=real64), allocatable :: freq(:)
 real(kind=real64), allocatable :: compfreq(:)
 real(kind=real64), allocatable :: compfreq_test(:)
 real(real64), allocatable :: B(:,:)
 real(real64), allocatable :: C(:,:)
 type(tgeno) :: mat
 type(tgeno) :: mat_ascii
 !$ real(kind=real64)::c_t1,c_t2,c_ttot
 !$ real(kind=real64)::t1,t2,ttot

 !C code
 integer(kind=int8), pointer, contiguous :: ivalt(:,:)
 integer(kind=int8), pointer, contiguous :: ivalt_transposed(:,:)
 integer(c_int) ::  c_snps, c_indiv
 integer(c_int) ::  c_ncol
 real(c_double), allocatable :: c_f(:)
 real(c_double), allocatable :: c_freq(:)
 real(c_double), allocatable :: c_B(:,:)
 real(c_double), allocatable :: c_C(:,:)
 type(c_ptr) :: c_plinkbed
 type(c_ptr) :: c_plinkbed_transposed
 type(c_ptr) :: c_compressed
 
 character(len=4) :: trans
 integer :: f_usegpu

 call get_command_argument(1,value=mat%genfile,status=io)
 call get_command_argument(2,value=freqfile,status=io)


#ifdef CUDA
 write(*,'("CUDA VERSION")')
#endif

 write(*,'(/2a)')' Bed file        : ',trim(mat%genfile)
 write(*,'(2a)')' Allele freq file: ',trim(freqfile)

 !READ FILES and ALLOCATE NEEDED ARRAYS
 call mat%read_bed()

 freq = readfreq(freqfile, mat%nsnp)

 !Input and output double precision arrays
 allocate(c_B(mat%nsnp, ncol))
 allocate(c_C(mat%ngen, ncol))
 call initarrays_c(c_B, c_C)

 if(lverbose)then
  do i = 1, 1
   write(*,'(a,*(1x,g0))')'f',c_B(i,1:4)
  enddo
 endif

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 allocate(B, source = real(c_B, real64))

 write(*,'("Leading dimension of B: ",i0)')size(B, 1)

 if(ltest)then
  mat_ascii%genfile = mat%genfile
  call mat_ascii%read_asci()

  allocate(compfreq_test(mat_ascii%nsnp))
  do i = 1, mat%nsnp
   compfreq_test(i) = sum(mat_ascii%cov%val(i,:), mat_ascii%cov%val(i,:).le.2._real64+epsilon(1._real64))&
                   /(count(mat_ascii%cov%val(i,:).le.2._real64+epsilon(1._real64)) * 2._real64)
  enddo

  !call check_freq('AR', compfreq_test, freq)

  !freq = compfreq_test

  deallocate(compfreq_test)

  if(lcenter)call center_by_f(mat_ascii, freq)
  allocate(C, mold = real(c_C, real64))
 endif

 !C code
 !INIT

#ifdef CUDA
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Transpose the 1-byte genotype matrix
 print*, mat%ngen, mat%nsnp, size(mat%cov%ival, 2),mat%cov%nrows
 call transpose_integermatrix(mat%ngen, mat%nsnp, size(mat%cov%ival, 2)&
                              , mat%cov%ival, idummy&
                              , trans_ival, compressedn, lda_atrans)
 write(*,'(a,i0,a,i0)')'Dimension of trans_ival: '&
                       , lda_atrans, ' x ', compressedn

 allocate(ivalt_transposed, source = transpose(trans_ival))
 c_plinkbed_transposed = c_loc(ivalt_transposed)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

 allocate(ivalt, source = transpose(mat%cov%ival))

 if(lverbose)then
  print*,'f byte 1:5,1 = ', ivalt(1:5,1)
  print'(a,*(1x, b8.8))','f byte 1:5,1 = ', ivalt(1:5,1)
 endif

 c_plinkbed = c_loc(ivalt)

 c_snps = int(mat%nsnp, c_int)
 c_indiv = int(mat%ngen, c_int)
 c_ncol = int(ncol, c_int)
 c_f = real(freq, c_double)

 !$ t1 = omp_get_wtime()
 call c_setOptions_compressed(usegpu, 0, 0, 0, 1, c_not_center, 0, 0, 256, c_lverbose)
 call c_plink2compressed(c_plinkbed, c_plinkbed_transposed, c_snps, c_indiv, c_f, c_ncol, c_compressed)
 !$ write(*,'(a,g0.6)')'  c_plink2compressed (s): ', omp_get_wtime() - t1

#ifndef CUDA
 if(ltest)then
  allocate(c_freq, mold = c_f)
  call c_get_freq(c_compressed, c_freq)
  call check_freq('CF', c_freq, freq)
  deallocate(c_freq)
 endif
#endif


 !MULTIPLICATION
 !Start repetition
 if(.not.lcenter)freq = 0._real64
 !$ c_ttot = 0
 !$ ttot = 0
 err = 0

 do irepet = 1, nrepet

  !C dgemm
  write(*,'(a)')"C dgemm"
!  trans = 'x'//c_null_char
  !$ c_t1 = omp_get_wtime()
  call c_dgemm_compressed('n', c_compressed, c_ncol, c_B, c_snps, c_C, c_indiv)
  !$ c_t2 = omp_get_wtime()
  !$ write(*,'(a,g0.6)')'  Elapsed time - C: ', c_t2 - c_t1
  !$ if (irepet > nmin) c_ttot = c_ttot + c_t2 - c_t1
  

  if(ltest)then

   !Compiler dgemm
   write(*,'(a)')"Compiler dgemm"
   !$ t1 = omp_get_wtime()
   C(:,:) = matmul(transpose(mat_ascii%cov%val), B(:,:))
   !$ t2 = omp_get_wtime()
   !$ write(*,'(a,g0.6)')'  Elapsed time - M: ', t2 - t1
   !$ if (irepet > nmin) ttot = ttot + t2 - t1

   !C check
   if (.true.) then 
     call check_matrix('M', 'C', C, c_C, 10, err)
   else
      write(*,'("No C testing !")')
   endif

  endif

 enddo
 !$ write(*,'("Average time - C: ", g0.6)')c_ttot/real(nrepet-nmin)
 !$ if(ltest)write(*,'("Average time - M: ", g0.6)')ttot/real(nrepet-nmin)

 if (err > 0) then  
   error stop 'Different outputs'
 endif

    !DEALLOCATE
 call c_free_compressed(c_compressed)

 deallocate(mat%cov%ival)
 deallocate(freq)
 deallocate(c_B, c_C)
 deallocate(B)

 if(ltest)then
  deallocate(C)
  deallocate(mat_ascii%cov%val)
 endif

contains

pure subroutine initarrays(B, C)
 real(real64), intent(out) :: B(:,:), C(:,:)

 integer :: i, j

 B = reshape([((-(10._real64*i+j), i=1, size(B, 1)), j=1, size(B, 2))]&
             , [size(B, 1), size(B, 2)])

 C = reshape([((-(10._real64*i+j), i=1, size(C, 1)), j=1, size(C, 2))]&
             , [size(C, 1), size(C, 2)])

end subroutine

pure subroutine initarrays_c(B, C)
 real(c_double), intent(out) :: B(:,:), C(:,:)
 
 integer :: i, j

 B = reshape([((-(10._c_double*i+j), i=1, size(B, 1)), j=1, size(B, 2))]&
             , [size(B, 1), size(B, 2)])

 C = reshape([((-(10._c_double*i+j), i=1, size(C, 1)), j=1, size(C, 2))]&
             , [size(C, 1), size(C, 2)])

end subroutine

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

subroutine check_freq(app, f, f_)
 character(*), intent(in) :: app
 real(real64), intent(in) :: f(:), f_(:)

 integer :: i

 if(lverbose)then
  do i =1, size(f)
   write(456,'(i0,1x,f10.8)')i,f(i)
  enddo
 endif

 if(any(abs(f_ - f) > 1.d-3))then
  do i = 1, size(f)
   if(abs(f_(i)-f(i)) > 1.d-3)then
    write(*, '(a,i0,*(1x,g0))')app//" freq ", i, f(i), f_(i)
   endif
  enddo
 else
  write(*,'(a)')app//': Allele frequencies: all abs differences < 10.^-3'
  write(*,'(a, g0)')app//': Max absolue difference: ',maxval(abs(f_ - f))
 endif

end subroutine

subroutine check_matrix(app_c, app_c_, C, C_, count_, err)
 character(*), intent(in) :: app_c, app_c_
 integer, intent(in) :: count_
 real(real64), intent(in) :: C(:,:)  !true
 real(real64), intent(in) :: C_(:,:) !to be checked
 integer, intent(inout) :: err

 integer :: i, j
 integer :: nrow, ncol
 integer :: counter
 character(:), allocatable :: app

 nrow = size(C, 1)
 ncol = size(C, 2)

 app = app_c//app_c_

 counter = 0

 if(any(abs(C - C_) > tol_real64))then      
  do j = 1, nrow
   do i = 1, ncol
    if (abs(C(j,i) - C_(j, i)) > 1.d-5) then
     if (counter < count_) then
      write(*,'(a, 2(1x, i0), *(1x, g0))')app, j, i, C(j,i), C_(j,i), C(j, i) - C_(j,i)
      counter = counter + 1
     endif
    endif
   enddo
  enddo
  write(*,'(a, g0.6)')app//' - Difference    : ', sum(C - C_)  
  write(*,'(a, g0.6)')app//' - abs Difference: ', sum(abs(C - C_))
  write(*,'(a, g0.6)')app//' - abs Difference of 1st row: ', sum(abs(C(1,:) - C_(1,:)))
  write(*,'(a, g0.6)')app//' - Difference of 1st col    : ', sum((C(:,1) - C_(:,1)))
  write(*,'(a, g0.6)')app//' - sum ('//app_c //')        : ', sum(C)  
  write(*,'(a, g0.6)')app//' - sum ('//app_c_//')        : ', sum(C_)  
  write(*,'(a, g0.6)')app//' - sum 1st col ('//app_c //'): ', sum(C(:,1))  
  write(*,'(a, g0.6)')app//' - sum 1st col ('//app_c_ //'): ', sum(C_(:,1))  
  err = 1
 else
  write(*,'("No '//app//' error !")')
 endif

end subroutine

end program
