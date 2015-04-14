module sparse_utils
! ==============================================================================
! Module for sparse matrix operations. It is partially a translation of F77
! code supplied with UMFPACK4.1. This code is written for the F-compiler, see
! http://www.fortran.com for the site. This compiler is available for free.
! As F is a subset of Fortran95, this code should run under any compliant
! Fortran95 compiler. Translation by Michael Janeschitz-Kriegl. Please report
! bugs to: info at heattransferconsult dot nl.  Please retain this header.
! ==============================================================================
! Although this code has been tested, it still may contain errors.
! Use with diligence! No warranty whatsoever. Use at own risk!
! ==============================================================================
use f90_kind

implicit none

private
public  :: triplet_read,triplet_write,triplet_spy,triplet_size, triplet_sort_col_cmplx
public  :: triplet_sort_col,triplet_sort_row,triplet_to_csc1,triplet_to_csc0
public  :: rua_read,rua_write,csc1_vector_mult,csc0_spy,csc1_spy
public  :: triplet_vector_mult,HBsize_read
public  :: csc_to_triplet_n0,count_non0
private :: triplet_swap, triplet_swap_cmplx

contains


!-------------------------------------------------------------------------------
!------------------------------------ READ -------------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_read(iunit,Ait,Ajt,Ax,ipoint)
  real(kind=double),dimension(:),intent(out)  :: Ax
  integer,dimension(:),intent(out)            :: Ait,Ajt,ipoint
  integer,intent(in) :: iunit
  integer  :: nz,i,ios

  nz = size(Ait)
  do i=1,nz
    ipoint(i) = i
    read(unit=iunit,fmt=*,iostat =ios) Ait(i),Ajt(i),Ax(i)
    if (ios /= 0) then
      print *,"Read error SR ""triplet_read"" at entry",i," of ",nz
      exit
    endif
  enddo

  end subroutine triplet_read


!-------------------------------------------------------------------------------
!------------------------------------ WRITE ------------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_write(iunit,Ait,Ajt,Ax,ipoint)
  real(kind=double),dimension(:),intent(in)  :: Ax
  integer,intent(in)                         :: iunit
  integer,dimension(:),intent(in)            :: Ait,Ajt
  integer,dimension(:),intent(inout)         :: ipoint
  integer  :: nz,i,ios

  nz = size(Ait)
  do i=1,nz
    write(unit=iunit,fmt="(2i5,es25.15,i5)",iostat=ios) &
    Ait(i),Ajt(i),Ax(i),ipoint(i)
    if (ios /= 0) then
      print *,"Write error SR ""triplet_write"" at entry",i," of ",nz
      exit
    endif
  enddo

  end subroutine triplet_write


!-------------------------------------------------------------------------------
!------------------------------------ SPYTRIPLET -------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_spy(iwrite,Ait,Ajt)
  integer,dimension(:),intent(in)             :: Ait,Ajt
  integer,intent(in)                          :: iwrite
  character(len=1),dimension(:,:),allocatable :: elemnt
  character(len=1),dimension(:),allocatable   :: header
  integer :: istat,i,j,n,nz
  character(len=3) :: string

  nz   = size(Ait)
  call triplet_size(Ait,Ajt,n)

  if (n > 1000) then
    write(unit=iwrite,fmt=*)"Matrix size>1000 (",n,"), print only &
                             &1000 x 1000 elements."
    n = 1000
  endif

  allocate(elemnt(n,n),header(n),stat=istat)
  if (istat /= 0) then
    write(unit=iwrite,fmt="(a50,i4)") &
   "Triplet-SPY: allocation error elemnt or header(n):",n
    stop
  endif

  elemnt = " "
  do i=1,nz
    elemnt(Ait(i),Ajt(i)) = "X"
  enddo

  header  = " "
  do i=1,n
    if (modulo(i,100) == 0) then
      write(unit=header(i),fmt="(i1)") modulo(i/100,100)
    endif
  enddo
  write(unit=iwrite,fmt="(a6,1000a1)") "      ",(header(i),i=1,n)
  header  = " "
  do i=1,n
    if (modulo(i,10) == 0) then
      write(unit=header(i),fmt="(i1)") modulo(i/10,10)
    endif
  enddo
  write(unit=iwrite,fmt="(a6,1000a1)") "      ",(header(i),i=1,n)
  write(unit=iwrite,fmt="(a6,1000i1)") "      ",(modulo(i,10),i=1,n)
  do i=1,n
    write(unit=string,fmt="(i3)") i
    write(unit=iwrite,fmt=*) string,"  ",(elemnt(i,j),j=1,n)
  enddo

  deallocate(elemnt,header,stat=istat)

  end subroutine triplet_spy



!-------------------------------------------------------------------------------
!------------------------------------ SIZE -------------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_size(Ait,Ajt,nsize)
  integer,dimension(:),intent(in) :: Ait,Ajt
  integer,intent(out)             :: nsize
  integer                         :: nz,i

  nz    = size(Ait)
  nsize = 0
  do i=1,nz
    if (Ait(i) > nsize) then
      nsize = Ait(i)
    endif
    if (Ajt(i) > nsize) then
      nsize = Ajt(i)
    endif
  enddo

  end subroutine triplet_size



!-------------------------------------------------------------------------------
!---------------------------------- SORT TO COLUMN -----------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_sort_col(Ait,Ajt,Ax,ipoint)
  real(kind=double),dimension(:),intent(inout)  :: Ax
  integer,dimension(:),intent(inout)            :: Ait,Ajt
  integer,dimension(size(Ax)),intent(out)       :: ipoint
  integer :: nz,i,j,k,icol_f,irow_f,irow_l

  nz=size(Ax)
! sort triplet form into columns, Ax(Ait(i),Ajt(i)), such that Ajt(i) ascends
  do j=1,nz-1
    do i=nz-1,j,-1
      if (Ajt(i) > Ajt(i+1)) then
        call triplet_swap(i,i+1,Ait,Ajt,Ax,ipoint)
      endif
    enddo
  enddo

  icol_f = Ajt(1)
  irow_f = 1
  do k=1,nz
    if (k==nz)then
      irow_l = nz
      do j=irow_f,irow_l-1
        do i=irow_l-1,irow_f,-1
          if(Ait(i) > Ait(i+1)) then
            call triplet_swap(i,i+1,Ait,Ajt,Ax,ipoint)
          endif
        enddo
      enddo
    endif
    if (Ajt(k) > icol_f) then
      irow_l = k-1
      do j=irow_f,irow_l-1
        do i=irow_l-1,irow_f,-1
          if(Ait(i) > Ait(i+1)) then
            call triplet_swap(i,i+1,Ait,Ajt,Ax,ipoint)
          endif
        enddo
      enddo
      icol_f = Ajt(k)
      irow_f = irow_l + 1
    endif
  enddo

  end subroutine triplet_sort_col

!-------------------------------------------------------------------------------
!---------------------------------- SORT TO COLUMN COMPLEX -----------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_sort_col_cmplx(Ait,Ajt,Ax,ipoint)
  complex,dimension(:),intent(inout)  :: Ax
  integer,dimension(:),intent(inout)            :: Ait,Ajt
  integer,dimension(size(Ax)),intent(out)       :: ipoint
  integer :: nz,i,j,k,icol_f,irow_f,irow_l

  nz=size(Ax)
! sort triplet form into columns, Ax(Ait(i),Ajt(i)), such that Ajt(i) ascends
  do j=1,nz-1
    do i=nz-1,j,-1
!print *, nz, j, i
      if (Ajt(i) > Ajt(i+1)) then
        call triplet_swap_cmplx(i,i+1,Ait,Ajt,Ax,ipoint)
      endif
    enddo
  enddo
print *, "1"
  icol_f = Ajt(1)
  irow_f = 1
  do k=1,nz
!print *, k
    if (k==nz)then
      irow_l = nz
      do j=irow_f,irow_l-1
        do i=irow_l-1,irow_f,-1
          if(Ait(i) > Ait(i+1)) then
            call triplet_swap_cmplx(i,i+1,Ait,Ajt,Ax,ipoint)
          endif
        enddo
      enddo
    endif
    if (Ajt(k) > icol_f) then
      irow_l = k-1
      do j=irow_f,irow_l-1
        do i=irow_l-1,irow_f,-1
          if(Ait(i) > Ait(i+1)) then
            call triplet_swap_cmplx(i,i+1,Ait,Ajt,Ax,ipoint)
          endif
        enddo
      enddo
      icol_f = Ajt(k)
      irow_f = irow_l + 1
    endif
  enddo

  end subroutine triplet_sort_col_cmplx



!-------------------------------------------------------------------------------
!---------------------------------- SORT TO ROW --------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_sort_row(Ait,Ajt,Ax,ipoint)
  real(kind=double),dimension(:),intent(inout)  :: Ax
  integer,dimension(:),intent(inout)            :: Ait,Ajt
  integer,dimension(size(Ax)),intent(out)       :: ipoint
  integer :: nz,i,j,k,icol_f,icol_l,irow_f

  nz=size(Ax)
! sort triplet form into rows, Ax(Ait(i),Ajt(i)), such that Ait(i) ascends
  do j=1,nz-1
    do i=nz-1,j,-1
      if (Ait(i) > Ait(i+1)) then
        call triplet_swap(i,i+1,Ait,Ajt,Ax,ipoint)
      endif
    enddo
  enddo

  irow_f = Ait(1)
  icol_f = 1
  do k=1,nz
    if (k==nz)then
      icol_l = nz
      do j=icol_f,icol_l-1
        do i=icol_l-1,icol_f,-1
          if(Ajt(i) > Ajt(i+1)) then
            call triplet_swap(i,i+1,Ait,Ajt,Ax,ipoint)
          endif
        enddo
      enddo
    endif
    if (Ait(k) > irow_f) then
      icol_l = k-1
      do j=icol_f,icol_l-1
        do i=icol_l-1,icol_f,-1
          if(Ajt(i) > Ajt(i+1)) then
            call triplet_swap(i,i+1,Ait,Ajt,Ax,ipoint)
          endif
        enddo
      enddo
      irow_f = Ait(k)
      icol_f = icol_l + 1
    endif
  enddo


  end subroutine triplet_sort_row



!-------------------------------------------------------------------------------
!------------------------------------ SWAP -------------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_swap(k,l,Ait,Ajt,Ax,ipoint)
  real(kind=double),dimension(:),intent(inout)  :: Ax
  integer,dimension(:),intent(inout)            :: Ait,Ajt,ipoint
  integer,intent(in)          :: k,l
  integer                     :: itemp
  real(kind=double)           :: dtemp

  itemp     = ipoint(k)
  ipoint(k) = ipoint(l)
  ipoint(l) = itemp

  itemp     = Ait(k)
  Ait(k)    = Ait(l)
  Ait(l)    = itemp

  itemp     = Ajt(k)
  Ajt(k)    = Ajt(l)
  Ajt(l)    = itemp

  dtemp     = Ax(k)
  Ax(k)     = Ax(l)
  Ax(l)     = dtemp

  end subroutine triplet_swap

!-------------------------------------------------------------------------------
!------------------------------------ SWAP COMPLEX------------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_swap_cmplx(k,l,Ait,Ajt,Ax,ipoint)
  complex,dimension(:),intent(inout)  :: Ax
  integer,dimension(:),intent(inout)            :: Ait,Ajt,ipoint
  integer,intent(in)          :: k,l
  integer                     :: itemp
  complex           :: dtemp

  itemp     = ipoint(k)
  ipoint(k) = ipoint(l)
  ipoint(l) = itemp

  itemp     = Ait(k)
  Ait(k)    = Ait(l)
  Ait(l)    = itemp

  itemp     = Ajt(k)
  Ajt(k)    = Ajt(l)
  Ajt(l)    = itemp

  dtemp     = Ax(k)
  Ax(k)     = Ax(l)
  Ax(l)     = dtemp

  end subroutine triplet_swap_cmplx


!-------------------------------------------------------------------------------
!---------------------------- TRIPLET-VECTOR MULTIPLY --------------------------
!-------------------------------------------------------------------------------
  subroutine triplet_vector_mult(Ait,Ajt,A,x,b)
! multiply  A * x = b
  integer,dimension(:),intent(in)            :: Ait,Ajt
  real(kind=double),dimension(:),intent(in)  :: A,x
  real(kind=double),dimension(:),intent(out) :: b
  integer :: nz,i
  nz = size(Ait)
  b = 0.0e0
  do i=1,nz
    b(Ait(i)) = b(Ait(i)) + A(i) * x(Ajt(i))
  enddo

  end subroutine triplet_vector_mult



!-------------------------------------------------------------------------------
!-------------------- CONVERT TO Compressed Sparse Column ----------------------
!-------------------------------------------------------------------------------
  subroutine triplet_to_csc1(Ait,Ajt,Ap,Ai)
! BEFORE CALLING THIS ROUTINE, THE TRIPLET FORM MATRIX
! MUST BE ORDERED ACCORDING TO COLUMNS. (using triplet_sort_col, above)
  integer,dimension(:),intent(in)  :: Ait,Ajt
  integer,dimension(:),intent(out) :: Ap,Ai
  integer :: nz,icol,istart,i,j,iend
  
  nz = size(Ait)
  icol   = 1
  istart = 1
  do i=1,nz
    if (Ajt(i) /= icol .or. i == nz) then
      Ap(icol) = istart
!print *,"icol,Ap(icol)",icol,Ap(icol)
      icol = icol + 1
      iend = i - 1
      if (i == nz) then
        iend = nz
      endif
      do j=istart,iend
        Ai(j) = Ait(j)
!print *,"                j,Ai(j)",j,Ai(j)
      enddo
      istart = i
    endif
  enddo
  Ap(icol) = nz + 1
!print *,"icol,Ap(icol)",icol,Ap(icol)
  
  end subroutine triplet_to_csc1



!-------------------------------------------------------------------------------
!-------------------- CONVERT TO Compressed Sparse Column C-format -------------
!-------------------------------------------------------------------------------
  subroutine triplet_to_csc0(Ait,Ajt,Ap,Ai)
! THE TRIPLET FORM MATRIX MUST BE ORDERED ACCORDING TO COLUMNS
! BEFORE CALLING THIS SUBROUTINE. (using triplet_sort_col, above)
  integer,dimension(:),intent(in)  :: Ait,Ajt
  integer,dimension(:),intent(out) :: Ap,Ai
  integer :: nz,icol,istart,i,j,iend
  
  nz = size(Ait)
  icol   = 1
  istart = 1
  do i=1,nz
    if (Ajt(i) /= icol .or. i == nz) then
      Ap(icol) = istart - 1
!print *,"icol,Ap(icol)",icol-1,Ap(icol)
      icol = icol + 1
      iend = i - 1
      if (i == nz) then
        iend = nz
      endif
      do j=istart,iend
        Ai(j) = Ait(j) - 1
!print *,"                j,Ai(j)",j-1,Ai(j)
      enddo
      istart = i
    endif
  enddo
  Ap(icol) = nz
!print *,"icol,Ap(icol)",icol-1,Ap(icol)
  
  end subroutine triplet_to_csc0



!-------------------------------------------------------------------------------
!--------------------------- CSC-MATRIX-VECTOR MULTIPLY ------------------------
!-------------------------------------------------------------------------------
  subroutine csc1_vector_mult(Ap,Ai,A,x,b)
! multiply  A * x = b
  integer,dimension(:),intent(in)            :: Ai,Ap
  real(kind=double),dimension(:),intent(in)  :: A,x
  real(kind=double),dimension(:),intent(out) :: b
  real(kind=double)  :: Aij
  integer            :: n,i,p,j

  n = size(Ap)-1
! Initialize b
  b = 0.0e0
! b = A * x
  do j = 1,n
    do p = Ap(j), Ap(j+1) - 1
      i    = Ai(p)
      Aij  = A(p)
      b(i) = b(i) + Aij * x(j)
    enddo
  enddo
  end subroutine csc1_vector_mult



!-------------------------------------------------------------------------------
!--------------------------- CSC1-MATRIX-SPY -----------------------------------
!-------------------------------------------------------------------------------
  subroutine csc1_spy(iwrite,Ap,Ai)
  integer,dimension(:),intent(in)             :: Ai,Ap
  integer,intent(in)                          :: iwrite
  character(len=1),dimension(:,:),allocatable :: elemnt
  character(len=1),dimension(:),allocatable   :: header
  character(len=3)   :: string
  integer            :: n,i,p,j,istat

  n = size(Ap)-1
  if (n > 1000) then
    write(unit=iwrite,fmt=*)"Matrix size>1000 (",n,"), print only &
                             &1000 x 1000 elements."
    n = 1000
  endif

! Allocate elemnt-array
  Allocate(elemnt(n,n),header(n),stat=istat)
  if (istat /= 0) then
    write(unit=iwrite,fmt="(a50,i4)") &
   "CSC1-SPY:    allocation error elemnt or header(n):",n
    stop
  endif

  elemnt=" "
  do j = 1,n
    do p = Ap(j), Ap(j+1) - 1
       i = Ai(p)
       if (i <= n) then
          elemnt(i,j) = "X"
       endif
    enddo
  enddo

  header = " "
  do i=1,n
    if (modulo(i,100) == 0) then
      write(unit=header(i),fmt="(i1)") modulo(i/100,100)
    endif
  enddo
  write(unit=iwrite,fmt="(a6,1000a1)") "      ",(header(i),i=1,n)

  header = " "
  do i=1,n
    if (modulo(i,10) == 0) then
      write(unit=header(i),fmt="(i1)") modulo(i/10,10)
    endif
  enddo
  write(unit=iwrite,fmt="(a6,1000a1)") "      ",(header(i),i=1,n)

  write(unit=iwrite,fmt="(a6,1000i1)") "      ",(modulo(i,10),i=1,n)
  do i=1,n
    write(unit=string,fmt="(i3)") i
    write(unit=iwrite,fmt=*) string,"  ",(elemnt(i,j),j=1,n)
  enddo

  deallocate(elemnt,header,stat=istat)

  end subroutine csc1_spy



!-------------------------------------------------------------------------------
!--------------------------- CSC0-MATRIX-SPY -----------------------------------
!-------------------------------------------------------------------------------
  subroutine csc0_spy(iwrite,Ap,Ai)
  integer,dimension(:),intent(in)             :: Ai,Ap
  integer,intent(in)                          :: iwrite
  character(len=1),dimension(:,:),allocatable :: elemnt
  character(len=1),dimension(:),allocatable   :: header
  character(len=3)   :: string
  integer            :: n,i,p,j,istat

  n = size(Ap)-1
! Allocate elemnt-array
  if (n > 1000) then
    write(unit=iwrite,fmt=*)"Matrix size>1000 (",n,"), print only &
                             &1000 x 1000 elements."
    n = 1000
  endif
  Allocate(elemnt(n,n),header(n),stat=istat)
  if (istat /= 0) then
    write(unit=iwrite,fmt="(a50,i4)") &
   "CSC1-SPY:    allocation error elemnt or header(n):",n
    stop
  endif

  elemnt=" "
  do j = 1,n
    do p = Ap(j), Ap(j+1) - 1
      i             = Ai(p+1)
      if (i <= n-1) then
        elemnt(i+1,j) = "X"
      endif
    enddo
  enddo

  header = " "
  do i=1,n
    if (modulo(i,100) == 0) then
      write(unit=header(i),fmt="(i1)") modulo(i/100,100)
    endif
  enddo
  write(unit=iwrite,fmt="(a6,1000a1)") "      ",(header(i),i=1,n)

  header = " "
  do i=1,n
    if (modulo(i,10) == 0) then
      write(unit=header(i),fmt="(i1)") modulo(i/10,10)
    endif
  enddo
  write(unit=iwrite,fmt="(a6,1000a1)") "      ",(header(i),i=1,n)

  write(unit=iwrite,fmt="(a6,1000i1)") "      ",(modulo(i,10),i=1,n)
  do i=1,n
    write(unit=string,fmt="(i3)") i
    write(unit=iwrite,fmt=*) string,"  ",(elemnt(i,j),j=1,n)
  enddo

  deallocate(elemnt,header,stat=istat)

  end subroutine csc0_spy



!-------------------------------------------------------------------------------
!-------------------------------- READ RUA-FORMAT ------------------------------
!-------------------------------------------------------------------------------
  subroutine rua_read(iread,Ap,Ai,Ax,title)

  real(kind=double),dimension(:),intent(out) :: Ax
  integer,dimension(:),intent(out) :: Ap,Ai
  integer,intent(in) :: iread
  integer :: nz, ios, p,                             &
             totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
             ncol, nrow, nel
  character(len=72),intent(out) :: title
  character(len=30) :: key
  character(len=3)  :: type
  character(len=16) :: ptrfmt
  character(len=16) :: indfmt
  character(len=20) :: valfmt
  character(len=20) :: rhsfmt

  read (unit=iread,fmt="(a72,a8)",iostat=ios) title, key
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix line 1"
    stop
  endif
  read (unit=iread,fmt="(5i14)",iostat=ios) totcrd,ptrcrd,indcrd,valcrd,rhscrd
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix line 2"
    stop
  endif
  read (unit=iread,fmt="(a3,tr11,4i14)",iostat=ios) type, nrow, ncol, nz, nel
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix line 3"
    stop
  endif
  read (unit=iread,fmt="(2a16,2a20)",iostat=ios) ptrfmt, indfmt, valfmt, rhsfmt
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix line 4"
    stop
  endif
  if (rhscrd > 0) then
      print *, "WARNING RHS of Harwell/Boeing matrix not read"
!      stop
  endif

  if (type /= "RUA" .or. nrow /= ncol) then
    print *,"RUA_READ: Error   can only handle square RUA matrices."
    stop
  endif

! read the matrix (0 or 1-based is determined by contents of Ap and Ai)
! thus, this routine can be used universally.
  read (unit=iread, fmt=ptrfmt,iostat=ios) (Ap (p), p = 1, ncol+1)
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix during Ap(p)"
    stop
  endif
  read (unit=iread, fmt=indfmt,iostat=ios) (Ai (p), p = 1, nz)
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix during Ai(p)"
    stop
  endif
  read (unit=iread, fmt=valfmt,iostat=ios) (Ax (p), p = 1, nz)
  if (ios /= 0)then
    print *, "Read error: Harwell/Boeing matrix during Ax(p)"
    stop
  endif

  end subroutine rua_read



!-------------------------------------------------------------------------------
!---------------------------- READ HB-MATRIX SIZE ------------------------------
!-------------------------------------------------------------------------------
  subroutine HBsize_read(iread,nrow,ncol,nz,rhs)
! ----------------------------------------------------------------------
! readhb_size:
!       read a sparse matrix in the Harwell/Boeing format and output the
!       size of the matrix (# rows, # columns, and # of entries)
! ----------------------------------------------------------------------
  integer,intent(in)  :: iread
  integer,intent(out) :: nrow, ncol, nz
  integer             :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, nel, ios
  character(len=72)   :: title
  character(len=30)   :: key
  character(len=16)   :: ptrfmt, indfmt
  character(len=20)   :: valfmt, rhsfmt
  character(len= 3)   :: type
  logical,intent(out) :: rhs

!-----------------------------------------------------------------------

! read header information from Harwell/Boeing matrix
  read(unit=iread,fmt="(a72, a8)",iostat=ios) title, key
  if (ios /= 0) then
    print *,"HBsize_read: read error in HB matrix header line 1."
    stop
  endif
  read(unit=iread,fmt="(5i14)"   ,iostat=ios) totcrd,ptrcrd,indcrd,valcrd,rhscrd
  if (ios /= 0) then
    print *,"HBsize_read: read error in HB matrix header line 2."
    stop
  endif
  read(unit=iread,fmt="(a3,tr11,4i14)",iostat=ios) type, nrow, ncol, nz, nel
  if (ios /= 0) then
    print *,"HBsize_read: read error in HB matrix header line 3."
    stop
  endif
  read(unit=iread,fmt="(2a16, 2a20)"  ,iostat=ios) ptrfmt, indfmt, valfmt, rhsfmt
  if (ios /= 0) then
    print *,"HBsize_read: read error in HB matrix header line 4."
    stop
  endif

  rhs = .false.
  if (rhscrd > 0) then
    rhs = .true.
  endif

  rewind(unit=iread)

  end subroutine HBsize_read



!-------------------------------------------------------------------------------
!-------------------------------- WRITE RUA-FORMAT ------------------------------
!-------------------------------------------------------------------------------
  subroutine rua_write(iwrite,Ap,Ai,Ax,title,key)

  real(kind=double),dimension(:),intent(in) :: Ax
  integer,dimension(:),intent(in) :: Ap,Ai
  integer,intent(in) :: iwrite
  integer :: n, nz,                                  &
             totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
             ncol, nrow, nel, ios, p
  character(len=72),intent(in) :: title
  character(len=8) ,intent(in) :: key
  character(len=3)  :: type   = "RUA"
  character(len=16) :: ptrfmt = "(16i5)          "
  character(len=16) :: indfmt = "(16i5)          "
  character(len=20) :: valfmt = "(3es25.15)          "
  character(len=20) :: rhsfmt = "                    "

  nz = size(Ax)
  n  = size(Ap)-1
  nrow = n
  ncol = n

! card1
  write (unit=iwrite,fmt="(a72,a8)",iostat=ios) title, key
  if (ios /= 0)then
     print *, "Write error: Harwell/Boeing matrix line 1"
     stop
  endif

! card2
  ptrcrd = (n+1)/16
  if (ptrcrd * 16 < n+1) then
    ptrcrd = ptrcrd + 1
  endif
  indcrd = nz/16 
  if (indcrd * 16 < nz) then
    indcrd = indcrd + 1
  endif
  valcrd = nz/3
  if (valcrd * 3 < nz) then
    valcrd = valcrd + 1
  endif
  totcrd = ptrcrd + indcrd + valcrd
  rhscrd = 0
  write (unit=iwrite,fmt="(5i14)",iostat=ios) totcrd,ptrcrd,indcrd,valcrd,rhscrd
  if (ios /= 0)then
     print *, "Write error: Harwell/Boeing matrix line 2"
     stop
  endif

! card3
  nel = 0
  write (unit=iwrite,fmt="(a3,tr11,4i14)",iostat=ios) type,nrow,ncol,nz,nel
  if (ios /= 0)then
     print *, "Write error: Harwell/Boeing matrix line 3"
     stop
  endif

! card4
  write (unit=iwrite,fmt="(2a16,2a20)",iostat=ios) ptrfmt,indfmt,valfmt,rhsfmt
  if (ios /= 0) then
     print *, "Write error: Harwell/Boeing matrix line 4"
     stop
  endif

! write the matrix (1 or 0-based)
  write (unit=iwrite, fmt=ptrfmt,iostat=ios) (Ap (p), p = 1, ncol+1)
  if (ios /= 0)then
     print *, "Write error: Harwell/Boeing matrix during Ap(p)"
     stop
  endif
  write (unit=iwrite, fmt=indfmt,iostat=ios) (Ai (p), p = 1, nz)
  if (ios /= 0)then
     print *, "Write error: Harwell/Boeing matrix during Ai(p)"
     stop
  endif
  write (unit=iwrite, fmt=valfmt,iostat=ios) (Ax (p), p = 1, nz)
  if (ios /= 0)then
     print *, "Write error: Harwell/Boeing matrix during Ax(p)"
     stop
  endif

  end subroutine rua_write



!-----------------------------------------------------------------------
!------------------ CONVERT CSC TO TRIPLET W/O ZEROES ------------------
!-----------------------------------------------------------------------
  subroutine csc_to_triplet_n0(Ap,Ai,Ax,Ait,Ajt,Axt)
  integer,dimension(:),intent(in)  :: Ap,Ai
  integer,dimension(:),intent(out) :: Ait,Ajt
  real(kind=double),dimension(:),intent(in)  :: Ax
  real(kind=double),dimension(:),intent(out) :: Axt
  integer :: n,nz,nnz,ne,p,icol,irow

! create the triplet form of the input matrix
! this works both for 0_based and 1_based arrays

  nz  = size(Ai)
  n   = size(Ap)-1
  call count_non0(Ax,nnz)

  ne  = 0
  do icol = 1, n
    do p  = Ap(icol), Ap(icol+1) - 1
      irow = Ai(p)
      if (Ax(p) /= 0.0) then
        ne = ne + 1
        Ait(ne) = irow
        Ajt(ne) = icol
        Axt(ne) = Ax(p)
      endif
    enddo
  enddo

  end subroutine csc_to_triplet_n0



!-----------------------------------------------------------------------
!------------------ COUNT NONZERO ENTRIES IN A VECTOR ------------------
!-----------------------------------------------------------------------
  subroutine count_non0(vector,ne)
  real(kind=double),dimension(:),intent(in) :: vector
  integer,intent(out) :: ne
  integer :: nz,i
  nz = size(vector)
  ne = 0
  do i=1,nz
    if(vector(i) /= 0.0)then
      ne=ne+1
    endif
  enddo
  end subroutine count_non0


!-------------------------------------------------------------------------------
!------------------------------------ END --------------------------------------
!-------------------------------------------------------------------------------

end module sparse_utils
