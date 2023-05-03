module modtestplink
 use, intrinsic:: iso_fortran_env, only: int8, int32, real64
 use modplink, only: plinkbedr
 implicit none

 private
 public::tgeno

 type::covariable
  integer(kind=int32)::nrows,ncol               !Number of rows and columns of val(:,:)
  integer(kind=int32)::neffectivelevel          !Number of effective levels (i.e., levels really present)
  integer(kind=int32),allocatable::indice(:)    !contains the column of the covariable at the position of its level ; dim=nlevel
  integer(kind=int8),pointer,contiguous::ival(:,:)    !compressed matrix containing coefficients per level; dim=ncoef(compressed)*(#level_really_present)
  real(kind=real64),pointer,contiguous::val(:,:)       !matrix containing coefficients per level; dim=ncoef*(#level_really_present)
  contains
  procedure::allocreal
  procedure::allocint
  final::dealloccov_scal,dealloccov_vect
 end type

 type tgeno
!   private
   integer(int32)::ngen,nsnp
   character(len=30)::genfile
   type(covariable)::cov
  contains
   private
   procedure, public ::read_asci=>readcovariable_plinkbed_to_asci
   procedure, public ::read_bed=>readcovariable_plinkbed_to_bed
 end type


contains

subroutine readcovariable_plinkbed_to_bed(geno)
 class(tgeno),intent(inout)::geno

 integer::indx,i,j,k,io,ncol,nlines,nbytes

 indx=index(geno%genfile,'.',back=.true.)
   !possibility that the file is a Plink file
   !Check if Plink bed file and associated file
!   inquire(file=geno%genfile,exist=lexist)
!   if(.not.lexist)then
!    write(*,'(a)')' ERROR: The following file does not exist:',adjustl(geno%genfile(:len_trim(geno%genfile)))
!    stop
!   endif
!   inquire(file=geno%genfile(1:indx-1)//'.bim',exist=lexist)
!   if(.not.lexist)then
!    write(*,'(a)')' ERROR: The following file does not exist:',adjustl(geno%genfile(1:indx-1)//'.bim')
!    stop
!   endif
!   inquire(file=geno%genfile(1:indx-1)//'.fam',exist=lexist)
!   if(.not.lexist)then
!    write(*,'(a)')' ERROR: The following file does not exist:',adjustl(geno%genfile(1:indx-1)//'.fam')
!    stop
!   endif
   !read .fam file for number of genotyped individuals
   open(newunit=i,file=geno%genfile(1:indx-1)//'.fam',action='read',status='old')!,buffered='yes')
   nlines=0
   do
    read(i,*,iostat=io)
    if(io.ne.0)exit
    nlines=nlines+1
   enddo
   close(i)
   !read .bim file for number of SNPs
   open(newunit=i,file=geno%genfile(1:indx-1)//'.bim',action='read',status='old')!,buffered='yes')
   ncol=0
   do
    read(i,*,iostat=io)
    if(io.ne.0)exit
    ncol=ncol+1
   enddo
   close(i)

   write(*,'(/"  Coefficient file: ",a)')adjustl(geno%genfile)

   nbytes=int(nlines/4)
   if(mod(nlines,4).ne.0)nbytes=nbytes+1


   !for efficiency, bed data are stored as opposite of other types of covariates
   allocate(geno%cov%ival(ncol,nbytes))
   geno%cov%ival=0
   geno%cov%neffectivelevel=nlines
   geno%cov%nrows=ncol                !number of SNPs
   geno%cov%ncol=nbytes            !number of animals

   !read .bed file
   open(newunit=i,file=geno%genfile,access='stream',status='old',action='read')!,buffered='yes')
   block
   integer(kind=int8)::tmp(3)
   read(i)tmp
   if(tmp(1).ne.108.or.tmp(2).ne.27.or.tmp(3).ne.1)then
    write(*,'(2a)')' ERROR: The following file is not a Plink SNP-major bed file: ',trim(geno%genfile)
    stop
   endif
   end block
   write(*,'(a,i0,a)')'   Number of SNPs: ', ncol
   write(*,'(a,i0,a)')'   Each SNP for all individuals requires ',nbytes,' bytes.'
   block
   integer::l,r,nread
   integer,parameter::READ_UNROLL=16
   integer(kind=int8)::tmp(nbytes,READ_UNROLL)
   do j=1,geno%cov%nrows,READ_UNROLL
     nread = min(READ_UNROLL, geno%cov%nrows-j+1)
     read(i,iostat=io)tmp(:,1:nread)
     if(io.ne.0)then
       write(*,'(2a)')'  ERROR unexpected while reading ', geno%genfile
       stop
     endif
     do l=1,nbytes
       !DIR$ UNROLL(READ_UNROLL)
       !DIR$ LOOP COUNT AVG(READ_UNROLL), MAX(READ_UNROLL)
       do r=1,nread
         geno%cov%ival(j+r-1,l) = tmp(l,r)
       enddo
     enddo
   enddo
   end block
   close(i)
 
   !for compatibility
   geno%nsnp = geno%cov%nrows
   geno%ngen = geno%cov%neffectivelevel

end subroutine

subroutine readcovariable_plinkbed_to_asci(geno)
 class(tgeno),intent(inout)::geno

 geno%ngen=-1
 geno%nsnp=-1
 call plinkbedr(geno%genfile,geno%cov%val,geno%ngen,geno%nsnp)

 write(*,'(/"  Number of genotypes :",i0)')geno%ngen
 write(*,'("  Number of snps      :",i0)')geno%nsnp

 geno%cov%neffectivelevel=geno%ngen
 geno%cov%nrows=geno%nsnp
 geno%cov%ncol=geno%ngen

end subroutine

subroutine allocreal(x,nrows,ncol)
 class(covariable),intent(inout)::x
 integer(kind=int32),intent(in)::nrows,ncol

 integer::stat

! if(allocated(x%val))deallocate(x%val)
 if(associated(x%val))then
  deallocate(x%val,stat=stat)
  if(stat.ne.0)nullify(x%val)
 endif
 x%nrows=nrows
 x%ncol=ncol
 allocate(x%val(x%nrows,x%ncol))
 x%val=0.d0

end subroutine

subroutine allocint(x,nrows,ncol)
 class(covariable),intent(inout)::x
 integer(kind=int32),intent(in)::nrows,ncol

 integer::stat

! if(allocated(x%val))deallocate(x%val)
 if(associated(x%val))then
  deallocate(x%val,stat=stat)
  if(stat.ne.0)nullify(x%val)
 endif
 x%nrows=nrows
 x%ncol=ncol
 allocate(x%ival(x%nrows,x%ncol))
 x%ival=0

end subroutine

!FINAL
subroutine dealloccov_scal(x)
 type(covariable)::x

 !integer::stat

 x%nrows=0
 x%ncol=0
 x%neffectivelevel=0

 !if(allocated(x%val))deallocate(x%val)
 if(associated(x%val))then
  !deallocate(x%val,stat=stat)
  !if(stat.ne.0)nullify(x%val)
  nullify(x%val)
 endif
 if(allocated(x%indice))deallocate(x%indice)
end subroutine

subroutine dealloccov_vect(x)
 type(covariable)::x(:)
 integer(kind=int32)::i

 !integer::stat

 x%nrows=0
 x%ncol=0
 x%neffectivelevel=0

 do i=1,size(x)
  !if(allocated(x(i)%val))deallocate(x(i)%val)
  if(associated(x(i)%val))then
   !deallocate(x(i)%val,stat=stat)
   !if(stat.ne.0)nullify(x(i)%val)
   nullify(x(i)%val)
  endif
  if(allocated(x(i)%indice))deallocate(x(i)%indice)
 enddo
end subroutine


end module
