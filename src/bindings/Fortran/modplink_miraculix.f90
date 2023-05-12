! License: MiXBLUP
!
! Authors 
! Jeremie Vandenplas, jeremie.vandenplas@wur.nl
!
! Copyright (C) 2022-2023 Jermie Vandenplas
!

module modplink_miraculix
 use, intrinsic:: iso_fortran_env, only:int8, int32, real32, real64
 !$ use omp_lib
 implicit none!(type, external)
 private
 public::plinkbedr
 public::transpose_integermatrix

 integer,parameter::wp=real64   !current precision for internal real computations

 !geno0 ! 00 = (0) = 0
 !geno1 ! 01 = (1) = 3     !missing value
 !geno2 ! 10 = (2) = 1
 !geno3 ! 11 = (3) = 2

 real(kind=wp),parameter::geno0=0._wp,geno1=3._wp,geno2=1._wp,geno3=2._wp   !geno1=missing
 real(kind=wp), parameter, public :: pgeno_wp(0:3) = [geno0, geno1, geno2, geno3]


 interface plinkbedr
  module procedure plinkbedr_r8
 end interface

 interface value2char
  module procedure inttochar_int32
  module procedure realtochar_real32
  module procedure realtochar_real64
 end interface

contains

function numidfam(filename) result(n)
 character(len=*), intent(in) :: filename

 integer :: n
 integer :: un

 n = 0
 open(newunit=un, file=trim(filename)//'.fam', status='old', action='read')
 call numlines(un,n)
 close(un)

end function

function numsnpbim(filename) result(n)
 character(len=*), intent(in) :: filename

 integer :: n
 integer :: un

 n = 0
 open(newunit=un, file=trim(filename)//'.bim', status='old', action='read')
 call numlines(un,n)
 close(un)

end function

subroutine plinkbedr_r8(namefile,geno,nids,nsnps,lfreq)
 integer(kind=int32),intent(inout),optional::nids,nsnps
 real(kind=real64),intent(out),pointer::geno(:,:)
 character(len=*),intent(in)::namefile
 logical,intent(in),optional::lfreq

 integer(kind=1)::tmp(3)
 integer(kind=int32)::i,j,io,un,nbytes,nan,nsnp
 real(kind=real64)::freq,tmpr
 !$ real(kind=real64)::t1
 character(len=:), allocatable :: filename
 logical::lfreq1

 i=index(namefile,'.',back=.true.)
 filename=namefile(1:i-1)
 write(*,'(a)')'Name of the file to append (without extension): '//trim(filename)

 if(.not.present(nids).or.(present(nids).and.nids.le.0))then
  nan=numidfam(filename)
 else
  nan=nids
 endif

 if(.not.present(nsnps).or.(present(nsnps).and.nsnps.le.0))then
  nsnp=numsnpbim(filename)
 else
  nsnp=nsnps
 endif
 write(*,'(a)')&
   'Number of animals and SNPs: '//value2char(nan)//' / '//value2char(nsnp)

 !if(allocated(geno))deallocate(geno)
 allocate(geno(nsnp,nan))

 !$ t1=omp_get_wtime()
 open(newunit=un,file=trim(filename)//'.bed',status='old',action='read',access='stream')
 read(un)tmp(:)
 if(tmp(1).eq.108.and.tmp(2).eq.27.and.tmp(3).eq.1)then
  nbytes=int(nan/4)
  if(mod(nan,4).ne.0)nbytes=nbytes+1

  block
  integer(1)::pos
  integer::jj,l,nread
!  integer,parameter::READ_UNROLL=1
  integer,parameter::READ_UNROLL=16
  integer(kind=1)::ttmp(nbytes,READ_UNROLL)
  real(kind=real64) :: geno_tmp(READ_UNROLL,4*nbytes)

  do j=1,nsnp,READ_UNROLL
    nread = min(READ_UNROLL, nsnp-j+1)
    read(un,iostat=io)ttmp(:,1:nread)
    if(io.ne.0)then
     write(*,'(a)')&
       'There is an issue during reading the Plink file.'// &
       'modplink '// &
       'plinkbedr_r8 ('//value2char(__LINE__)//')'
     error stop io
    endif
    do jj=1,nbytes
     !DIR$ LOOP COUNT AVG(READ_UNROLL), MAX(READ_UNROLL)
     do l=1,nread
      do pos = 0, 3
       geno_tmp(l, (jj-1)*4+1+pos) = pgeno_wp(ibits(ttmp(jj,l), pos*2 ,2))
      enddo
     enddo
    enddo
    geno(j:j+nread-1, 1:nan) = geno_tmp(1:nread, 1:nan)
  enddo
  end block

  close(un)
 else
  write(*,'(a)')&
    'It is not a SNP-major bed file of plink!'// &
    'modplink '// &
    'plinkbedr_r8 ('//value2char(__LINE__)//')'
  close(un)
  error stop
 endif
 !$ write(*,'(a,i0,a,f0.3,a)')'   Number of read SNP: ',i,' in ',omp_get_wtime()-t1,'s'

 !Recode for the minor allele if requested
 lfreq1=.false.
 if(present(lfreq))lfreq1=lfreq
 if(lfreq1)then
  !$ t1=omp_get_wtime()
  do i=1,nsnp
   freq=0.d0;io=0
   do j=1,nan
    tmpr=geno(i,j)
    if(tmpr<geno1)then
     freq=freq+tmpr
     io=io+1
    endif
   enddo
   freq=freq/(2*io)
   if(freq>=5.d-1)then
    do j=1,nan
     if(geno(i,j).le.geno3)geno(i,j)=geno3-geno(i,j)
    enddo
   endif
  enddo
  !$ write(*,'(a,f0.3,a)')'  Recoded SNP in: ',omp_get_wtime()-t1,'s'
 endif

 if(present(nids))nids=nan
 if(present(nsnps))nsnps=nsnp

end subroutine

subroutine transpose_integermatrix(m,k,k1,a,lda,atrans,compressedn,lda_atrans)
 integer(kind=int32),intent(in)::m,k,k1,lda
 integer(kind=int32),intent(out)::lda_atrans
 integer(kind=int8),intent(in)::a(:,:)
 integer(kind=int8),allocatable,intent(out)::atrans(:,:)

 integer(kind=int32)::i,j,l,compressedn,rest_snp
 integer(kind=int32)::pos_an,pos_snp,pos_snp_convert,min_pos_an
 real(kind=real64)::t1

 !$ t1=omp_get_wtime()
 compressedn=int(k/4)
 if(mod(k,4).ne.0)compressedn=compressedn+1
 allocate(atrans(m,compressedn))
 lda_atrans=m
 atrans=0
!print*,'bbbb a',size(a,1),size(a,2),' atrans ',size(atrans,1),size(atrans,2)

 do i=1,k1              !compressed animal
  min_pos_an=(i-1)*4
  do j=1,k              !snp
   pos_snp=j
   do l=1,4             !genotype in compressed format
    !retrieve position of animal and of snp
    pos_an=min_pos_an+l
    if(pos_an>m)cycle
    pos_snp_convert=int(pos_snp/4)
    rest_snp=4
    if(mod(pos_snp,4).ne.0)then
     pos_snp_convert=pos_snp_convert+1
     rest_snp=mod(pos_snp,4)
    endif
!print*,'zz',rest_snp,pos_snp/4
!print*,'xx',j,i,(l-1)*2
!print*,'yy',pos_an,pos_snp_convert,(rest_snp-1)*2
    call mvbits(a(j,i),(l-1)*2,2,atrans(pos_an,pos_snp_convert),(rest_snp-1)*2)
   enddo
  enddo
 enddo

 !$ write(*, '(a,f0.3)')'      Wall clock time (s) for transposition  :',omp_get_wtime()-t1
end subroutine



!!!!
subroutine numlines(unfile,n)
 integer(kind=int32)::io
 integer(kind=int32),intent(in)::unfile
 integer(kind=int32),intent(out)::n
 character(len=20)::a
 rewind(unfile)
 n=0
 do
  read(unfile,*,iostat=io)a
  if (io.ne.0) exit
  n=n+1
 enddo
 rewind(unfile)
end subroutine

subroutine numcol(unfile,n)
 integer(kind=int32),intent(in)::unfile
 character(len=1000000)::a
 integer(kind=int32),intent(out)::n
 integer(kind=int32)::curr,first,last,lena,i
 rewind(unfile)
 read(unfile,"(a)")a
 curr=1;lena=len(a);n=0
 do
  first=0
  do i=curr,lena
   if (a(i:i) /= " ") then
    first=i
    exit
   endif
  enddo
  if (first == 0) exit
  curr=first+1
  last=0
  do i=curr,lena
   if (a(i:i) == " ") then
    last=i
    exit
   endif
  enddo
  if (last == 0) last=lena
  n=n+1
  curr=last+1
 enddo
 rewind(unfile)
end subroutine

pure function inttochar_int32(i) result(res)
 character(:), allocatable :: res
 integer(int32), intent(in) :: i

 character(range(i)+2) :: tmp

 write(tmp,'(i0)') i
 res = trim(tmp)

end function

pure function realtochar_real32(i,fmt) result(res)
 character(:), allocatable :: res
 real(real32), intent(in) :: i

 character(*),intent(in), optional :: fmt

 character(range(i)+2) :: tmp
 character(len=:), allocatable :: fmt_

 fmt_ = '(g0)'
 if(present(fmt)) fmt_ = fmt

 write(tmp, fmt_) i

 res = trim(tmp)

end function
pure function realtochar_real64(i,fmt) result(res)
 character(:), allocatable :: res
 real(real64), intent(in) :: i

 character(*),intent(in), optional :: fmt

 character(range(i)+2) :: tmp
 character(len=:), allocatable :: fmt_

 fmt_ = '(g0)'
 if(present(fmt)) fmt_ = fmt

 write(tmp, fmt_) i

 res = trim(tmp)

end function

end module
