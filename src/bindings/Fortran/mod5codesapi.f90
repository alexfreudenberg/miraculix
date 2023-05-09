! License: MiXBLUP
!
! Authors 
! Jeremie Vandenplas, jeremie.vandenplas@wur.nl
!
! Copyright (C) 2022-2023 Jermie Vandenplas
!


module mod5codesapi
 use, intrinsic:: iso_c_binding, only: c_char, c_int, c_double, c_ptr
 implicit none
 private
 public :: c_setoptions_compressed
 public :: c_get_compressed_freq
 public :: c_plink2compressed
 public :: c_dgemm_compressed
 public :: c_free_compressed

 interface
  subroutine c_setoptions_compressed(use_gpu&
                                     , cores&
                                     , floatLoop&
                                     , meanSubstract&
                                     , ignore_missings&
                                     , do_not_center&
                                     , do_normalize&
                                     , use_miraculix_freq&
                                     , variant&
                                     , print_details)&
             bind(C, name='setOptions_compressed')
   import c_int
   integer(c_int), value, intent(in) :: use_gpu
   integer(c_int), value, intent(in) :: cores
   integer(c_int), value, intent(in) :: floatLoop
   integer(c_int), value, intent(in) :: meanSubstract
   integer(c_int), value, intent(in) :: ignore_missings
   integer(c_int), value, intent(in) :: do_not_center
   integer(c_int), value, intent(in) :: do_normalize
   integer(c_int), value, intent(in) :: use_miraculix_freq
   integer(c_int), value, intent(in) :: variant
   integer(c_int), value, intent(in) :: print_details
  end subroutine

  subroutine c_plink2compressed(plink, plink_transposed, snps, indiv, f, n, compressed)&
             bind(C, name='plink2compressed')
   import c_int, c_ptr, c_double
   type(c_ptr), value, intent(in) :: plink
   type(c_ptr), value, intent(in) :: plink_transposed
   integer(c_int), value, intent(in) :: snps
   integer(c_int), value, intent(in) :: indiv
   real(c_double), intent(in) :: f(*)
   integer(c_int), value, intent(in) :: n
   type(c_ptr), intent(out) :: compressed !(matches void**)
   !type(c_ptr), value, intent(in) :: compressed !(matches void*)
  end subroutine
 
  subroutine c_dgemm_compressed(trans, compressed, n, B, ldb, C, ldc)&
             bind(C, name='dgemm_compressed')
   import c_char, c_int, c_double, c_ptr
   character(c_char), intent(in) :: trans(*)
   type(c_ptr), value, intent(in) :: compressed !(matches void*)
   integer(c_int), value, intent(in) :: n
   real(c_double), intent(in) :: B(ldb, *)
   integer(c_int), value, intent(in) :: ldb
   real(c_double), intent(inout) :: C(ldc, *)
   integer(c_int), value, intent(in) :: ldc
  end subroutine

  subroutine c_get_compressed_freq(compressed, freq)&
             bind(C, name='get_compressed_freq')
   import c_double, c_ptr
   type(c_ptr), value, intent(in) :: compressed !(matches void*)
   real(c_double), intent(out) :: freq(*)
  end subroutine

  subroutine c_free_compressed(compressed)&
             bind(C, name='free_compressed')
   import c_ptr
   type(c_ptr), intent(inout) :: compressed !(matches void**)
  end subroutine
 end interface

contains
end module
