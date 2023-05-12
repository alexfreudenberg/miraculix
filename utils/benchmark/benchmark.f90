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

! This file benchmarks the 5codes algorithm and the GPU implementation 
! against the  original code by Vandenplas et al (2020) and the Intel MKL
! dgemm function
! A CUDA installation and a GPU of compute capability greater than 7.0 is required

program benchmark
 use, intrinsic:: iso_fortran_env, only: int8, real64
 use, intrinsic:: iso_c_binding, only: c_loc, c_int, c_ptr, c_null_ptr, c_double&
                                       , c_char, c_null_char
 use mod5codesapi, only: c_plink2compressed, c_dgemm_compressed, c_free_compressed&
                         , c_get_freq => c_get_compressed_freq&
                         , c_setOptions_compressed
 use modtestplink, only: tgeno
 use modplink_miraculix, only: dgemm_compressed, transpose_integermatrix
 use modplink_miraculix, only: allelefreq
 !$ use omp_lib
 implicit none

 integer, parameter :: ncol =  10
 integer, parameter :: nrepet = 10
 integer, parameter :: nmin = merge(1, 0, nrepet > 1)
 logical, parameter :: ltest = .true.
 logical, parameter :: lcenter = .true.
 logical, parameter :: lverbose = .true.

 character(len=3) :: test_mode
 character(len=1) :: line_break = char(10)


 integer(c_int), parameter :: c_not_center = merge(0, 1, lcenter)
 integer(c_int), parameter :: c_lverbose = merge(1, 0, lverbose)
 
 integer :: irepet
 integer :: i, j, io, err, nrow

 integer :: arg_count
 integer :: usegpu
 integer :: idummy, compressedn, lda_atrans
 integer(kind=int8), allocatable :: trans_ival(:,:)
 character(len=60) :: freqfile
 real(kind=real64), allocatable :: freq(:)
 real(real64), allocatable :: val_transposed(:,:)
 real(real64), allocatable :: B(:,:)
 real(real64), allocatable :: C(:,:)
 real(real64), allocatable :: f_C(:,:)
 real(real64), allocatable :: B_T(:,:)
 real(real64), allocatable :: C_T(:,:)
 real(real64), allocatable :: f_C_T(:,:)
 type(tgeno) :: mat
 type(tgeno) :: mat_ascii
 !$ real(kind=real64)::cpu_t1,cpu_t2,cpu_ttot
 !$ real(kind=real64)::gpu_t1,gpu_t2,gpu_ttot
 !$ real(kind=real64)::f_t1,f_t2,f_ttot
 !$ real(kind=real64)::mkl_t1,mkl_t2,mkl_ttot
 !$ real(kind=real64)::t1,t2,ttot
 real :: start_time, end_time, elapsed_time


 !C code
 integer(kind=int8), pointer, contiguous :: ivalt(:,:)
 integer(kind=int8), pointer, contiguous :: ivalt_transposed(:,:)
 integer(c_int) ::  c_snps, c_indiv
 integer(c_int) ::  c_ncol
 real(c_double), allocatable :: c_f(:)
 real(c_double), allocatable :: c_freq(:)
 real(c_double), allocatable :: c_B(:,:)
 real(c_double), allocatable :: c_C(:,:)
 real(c_double), allocatable :: c_B_T(:,:)
 real(c_double), allocatable :: c_C_T(:,:)
 type(c_ptr) :: c_plinkbed
 type(c_ptr) :: c_plinkbed_transposed
 type(c_ptr) :: c_compressed
 
! Get input arguments
 arg_count = command_argument_count()
 if (arg_count .ne. 3) then
  print *, "Program requires three input arguments: test mode, bed file, freq file."
  print *, "Test mode allows options 'Org' (Fortran implementation), 'CPU', 'GPU' and 'MKL'."
  stop
 end if

 call get_command_argument(1,value=test_mode,status=io)
 call get_command_argument(2,value=mat%genfile,status=io)
 call get_command_argument(3,value=freqfile,status=io)
 
 write(*,'(/2a)')' Test mode       : ',trim(test_mode)
 write(*,'(2a)')' Bed file        : ',trim(mat%genfile)
 write(*,'(2a)')' Allele freq file: ',trim(freqfile)

 !READ FILES and ALLOCATE NEEDED ARRAYS
 call mat%read_bed()

 freq = readfreq(freqfile, mat%nsnp)

 !Input and output double precision arrays
 allocate(c_B(mat%nsnp, ncol))
 allocate(c_C(mat%ngen, ncol))
 call initarrays_c(c_B, c_C)

 allocate(c_B_T(mat%ngen, ncol))
 allocate(c_C_T(mat%nsnp, ncol))
 call initarrays_c(c_B_T, c_C_T)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 allocate(B, source = real(c_B, real64))
 allocate(B_T, source = real(c_B_T, real64))
 allocate(f_C, mold = real(c_C, real64))
 allocate(f_C_T, mold = real(c_C, real64))
 
 allocate(ivalt, source = transpose(mat%cov%ival))

 !C code
 !INIT



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 c_plinkbed = c_loc(ivalt)
 c_snps = int(mat%nsnp, c_int)
 c_indiv = int(mat%ngen, c_int)
 c_ncol = int(ncol, c_int)
 c_f = real(freq, c_double)


 !MULTIPLICATION
 !Start repetition
 !$ cpu_ttot = 0
 !$ f_ttot = 0
 !$ mkl_ttot = 0
 !$ ttot = 0


 select case (test_mode)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!F dgemm
 case ("Org")
  write(*,'(5A)') line_break, line_break, line_break, line_break, line_break
  write(*,'(a)')"F dgemm"
  do irepet = 1, (nrepet+ nmin)
    !$ f_t1 = omp_get_wtime()
    call dgemm_compressed('t', 'n', mat%cov%neffectivelevel, ncol, mat%cov%nrows, mat%cov%ncol&
                          ,1._real64, mat%cov%ival, freq, mat%cov%nrows, B, mat%nsnp&
                          ,0._real64, f_C, mat%ngen)
    !$ f_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z - F: ', f_t2 - f_t1
    !$ if (irepet > nmin) f_ttot = f_ttot + f_t2 - f_t1
  enddo
  !$ write(*,'("Average time - Z - F: ", g0.6)')f_ttot/real(nrepet)
  !$ f_ttot = 0

  do irepet = 1, (nrepet+ nmin)
    !$ f_t1 = omp_get_wtime()
    call dgemm_compressed('n', 'n', mat%cov%nrows, ncol, mat%cov%neffectivelevel, mat%cov%ncol&
                          ,1._real64, mat%cov%ival, freq, mat%cov%nrows, B_T, mat%ngen&
                          ,0._real64, f_C_T, mat%nsnp)
    !$ f_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z^T - F: ', f_t2 - f_t1
    !$ if (irepet > nmin) f_ttot = f_ttot + f_t2 - f_t1
  enddo
  !$ write(*,'("Average time - Z^T - F: ", g0.6)')f_ttot/real(nrepet)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CPU dgemm
 case ("CPU")
  !$ t1 = omp_get_wtime()
  usegpu = 0
  c_plinkbed_transposed = c_plinkbed
  call c_setOptions_compressed(usegpu, 0, 0, 0, 1, c_not_center, 0, 0, 512, 1)
  call c_plink2compressed(c_plinkbed, c_plinkbed_transposed, c_snps, c_indiv, c_f, c_ncol, c_compressed)
  !$ write(*,'(a,g0.6)')'  c_plink2compressed (s): ', omp_get_wtime() - t1
  
  write(*,'(5A)') line_break, line_break, line_break, line_break, line_break
  write(*,'(a)')"CPU dgemm"
  do irepet = 1, (nrepet+nmin)
    !$ cpu_t1 = omp_get_wtime()
    call c_dgemm_compressed('n', c_compressed, c_ncol, c_B, c_snps, c_C, c_indiv)
    !$ cpu_t2 = omp_get_wtime()
    !$ if (irepet > nmin)  write(*,'(a,g0.6)')'  Elapsed time - Z - CPU: ', cpu_t2 - cpu_t1
    !$ if (irepet > nmin) cpu_ttot = cpu_ttot + cpu_t2 - cpu_t1
  enddo
  !$ write(*,'("Average time - Z - CPU: ", g0.6)')cpu_ttot/real(nrepet)
  !$ cpu_ttot = 0

  do irepet = 1, (nrepet+nmin)
    !$ cpu_t1 = omp_get_wtime()
    call c_dgemm_compressed('t', c_compressed, c_ncol, c_B_T, c_indiv, c_C_T, c_snps)  
    !$ cpu_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z^T - CPU: ', cpu_t2 - cpu_t1
    !$ if (irepet > nmin) cpu_ttot = cpu_ttot + cpu_t2 - cpu_t1
  enddo
  !$ write(*,'("Average time - Z^T - CPU: ", g0.6)')cpu_ttot/real(nrepet)
  !DEALLOCATE
  call c_free_compressed(c_compressed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GPU dgemm
 case ("GPU")
  !$ t1 = omp_get_wtime()
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
  usegpu = 1
  call c_setOptions_compressed(usegpu, 0, 0, 0, 1, c_not_center, 0, 0, 32, 1)
  call c_plink2compressed(c_plinkbed, c_plinkbed_transposed, c_snps, c_indiv, c_f, c_ncol, c_compressed)
  !$ write(*,'(a,g0.6)')'  c_plink2compressed (s): ', omp_get_wtime() - t1

  write(*,'(5A)') line_break, line_break, line_break, line_break, line_break
  write(*,'(a)')"GPU dgemm"
  do irepet = 1, (nrepet+nmin)
    !$ gpu_t1 = omp_get_wtime()
    call c_dgemm_compressed('n', c_compressed, c_ncol, c_B, c_snps, c_C, c_indiv)
    !$ gpu_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z - GPU: ', gpu_t2 - gpu_t1
    !$ if (irepet > nmin) gpu_ttot = gpu_ttot + gpu_t2 - gpu_t1
  enddo
  !$ write(*,'("Average time - Z - GPU: ", g0.6)')gpu_ttot/real(nrepet)
  !$ gpu_ttot = 0

  do irepet = 1, (nrepet+nmin)
    !$ gpu_t1 = omp_get_wtime()
    call c_dgemm_compressed('t', c_compressed, c_ncol, c_B_T, c_indiv, c_C_T, c_snps)  
    !$ gpu_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z^T - GPU: ', gpu_t2 - gpu_t1
    !$ if (irepet > nmin) gpu_ttot = gpu_ttot + gpu_t2 - gpu_t1
  enddo
  !$ write(*,'("Average time - Z^T - GPU: ", g0.6)')gpu_ttot/real(nrepet)
  !DEALLOCATE
  call c_free_compressed(c_compressed)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MKL dgemm
 case ("MKL")
  mat_ascii%genfile = mat%genfile
  call mat_ascii%read_asci()
  if(lcenter)call center_by_f(mat_ascii, freq)

  allocate(val_transposed, source = transpose(mat_ascii%cov%val))
  allocate(C, mold = real(c_C, real64))
  allocate(C_T, mold = real(c_C_T, real64))

  write(*,'(5A)') line_break, line_break, line_break, line_break, line_break
  write(*,'(a)')"MKL dgemm"

  !$ mkl_ttot = 0
  do irepet = 1, (nrepet+nmin)
    !$ mkl_t1 = omp_get_wtime()    
    C(:,:) = matmul(val_transposed, B(:,:))
    !$ mkl_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z - MKL: ', mkl_t2 - mkl_t1
    !$ if (irepet > nmin) mkl_ttot = mkl_ttot + mkl_t2 - mkl_t1
  enddo
  !$ write(*,'("Average time - Z - MKL: ", g0.6)')mkl_ttot/real(nrepet)

  !$ mkl_ttot = 0
  do irepet = 1, (nrepet+nmin)
    !$ mkl_t1 = omp_get_wtime()    
    C_T(:,:) =  matmul(mat_ascii%cov%val, B_T(:,:))
    !$ mkl_t2 = omp_get_wtime()
    !$ if (irepet > nmin) write(*,'(a,g0.6)')'  Elapsed time - Z^T - MKL: ', mkl_t2 - mkl_t1
    !$ if (irepet > nmin) mkl_ttot = mkl_ttot + mkl_t2 - mkl_t1
  enddo
  !$ write(*,'("Average time - Z^T - MKL: ", g0.6)')mkl_ttot/real(nrepet)

  deallocate(C, C_T)
  deallocate(mat_ascii%cov%val)
  
 case default
   print *, "Invalid option provided."
 end select

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 deallocate(mat%cov%ival)
 deallocate(freq)
 deallocate(c_B, c_C)
 deallocate(B, B_T)
 deallocate(c_B_T, c_C_T)
 deallocate(f_C_T, f_C)


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


end program
