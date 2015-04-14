module my_sparse_aggregate

     use precisions, only: wp, cwp, sp, dp
     use my_sparse
!     use save_load

!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!---------------------------- CONVERT coo2fullL SPARSE ---------------------------------------
     interface coo2full
          module procedure conv_coo2full
          module procedure conv_coo2full_int
          module procedure conv_coo2full_cmplx
     end interface
!---------------------------- CONVERT COO2CSR SPARSE ---------------------------------------
     interface coo2csr
          module procedure conv_coo2csr
          module procedure conv_coo2csr_dp
          module procedure conv_coo2csr_cmplx
          module procedure conv_coo2csr_cmplx_cmplx
     end interface  
!---------------------------- CONVERT CSR2COO SPARSE ----------------------------------
     interface csr2coo
          module procedure conv_csr2coo
          module procedure conv_csr2coo_dp
          module procedure conv_csr2coo_cmplx
     end interface  
!---------------------------- CONVERT CSR2CSC SPARSE ----------------------------------
     interface csr2csc					! For UMFPACK solver, also converts to 0-based indexing
!          module procedure conv_csr2csc
          module procedure conv_csr2csc_cmplx
     end interface  

!---------------------------- MULT DIAG BY SPARSE -----------------------------------------
     interface dia_csr_mul					! Better to use skit method with this function
          module procedure mul_dia_csr
          module procedure mul_dia_csr_cmplx
          module procedure mul_dia_int_csr_cmplx
     end interface
!---------------------------- CSR PLUS CSR ---------------------------------
     interface csr_csr_add					! SparseKit method not implememnted! Requires careful pre-allocation, tedious...
          module procedure add_csr
          module procedure add_csr_cmplx
     end interface
!---------------------------- CSR PLUS DIAG ---------------------------------
     interface csr_dia_add
          module procedure add_diacsr_cmplx	! NO MEMORY ALLOCATION. SparseKit only
     end interface
!==========================================================================================

contains
!==========================================================================================
subroutine conv_coo2full(mat_coo, mat_ful)

implicit none

     type (triplet)			:: mat_coo
     real(wp),allocatable	:: mat_ful(:,:)
	 integer :: istat, c

    allocate(mat_ful(mat_coo%ni, mat_coo%nj), stat = istat)

    mat_ful = 0

    do c = 1,mat_coo%nz
    	mat_ful(mat_coo%indi(c),mat_coo%indj(c)) = mat_ful(mat_coo%indi(c),mat_coo%indj(c)) + mat_coo%vals(c)
    enddo

end subroutine conv_coo2full

!==========================================================================================
subroutine conv_coo2full_int(mat_coo, mat_ful)

implicit none

     type (triplet_int)			:: mat_coo
     integer, allocatable		:: mat_ful(:,:)
	 integer :: istat, c

    allocate(mat_ful(mat_coo%ni, mat_coo%nj), stat = istat)

    mat_ful = 0

    do c = 1,mat_coo%nz
    	mat_ful(mat_coo%indi(c),mat_coo%indj(c)) = mat_ful(mat_coo%indi(c),mat_coo%indj(c)) + mat_coo%vals(c)
    enddo

end subroutine conv_coo2full_int

!==========================================================================================
subroutine conv_coo2full_cmplx(mat_coo, mat_ful)

implicit none

     type (triplet_cmplx)		:: mat_coo
     complex(cwp),allocatable	:: mat_ful(:,:)
	 integer :: istat, c

    allocate(mat_ful(mat_coo%ni, mat_coo%nj), stat = istat)

    mat_ful = 0

    do c = 1,mat_coo%nz
    	mat_ful(mat_coo%indi(c),mat_coo%indj(c)) = mat_ful(mat_coo%indi(c),mat_coo%indj(c)) + mat_coo%vals(c)
    enddo

end subroutine conv_coo2full_cmplx
!==========================================================================================

!==========================================================================================
subroutine conv_coo2csr(mat_coo, mat_csr, method)

implicit none

     type (triplet) :: mat_coo
     type (csr) :: mat_csr
     character :: method

     integer,dimension(:),allocatable :: iwork
     integer :: info_mkl, job(8), istat

     call init_sparse_0(mat_csr, mat_coo%ni, mat_coo%nj, mat_coo%nz)

if (method == skit) then
!print *, "mat_coo", mat_coo%ni, mat_coo%nj, mat_coo%nz
!sparskit convert COO to CSR:
     call coocsr(mat_coo%ni,mat_coo%nz,mat_coo%vals,mat_coo%indi, mat_coo%indj, &
                       mat_csr%vals,mat_csr%indj,mat_csr%indi)
!print *, "mat_csr", mat_csr%ni, mat_csr%nj, mat_csr%nz
!sparskit sort CSR with column index incr within each row:
     allocate(iwork(2*mat_csr%nz),stat=istat)
     call csort(mat_csr%ni,mat_csr%vals,mat_csr%indj,mat_csr%indi,iwork,.true.)

else

! MKL convert triplet/COO to CSR format
job(1)=2 ! convert COO to CSR and sort
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo%nz ! # of non-zeros in matrix A
job(6)=0 ! Fill all CSR arrays

if (wp==dp) then
     call mkl_dcsrcoo(job, mat_coo%ni, mat_csr%vals, mat_csr%indj, mat_csr%indi, mat_csr%nz, &
                 mat_coo%vals, mat_coo%indi, mat_coo%indj, info_mkl)
elseif (wp==sp) then
     call mkl_scsrcoo(job, mat_coo%ni, mat_csr%vals, mat_csr%indj, mat_csr%indi, mat_csr%nz, &
                 mat_coo%vals, mat_coo%indi, mat_coo%indj, info_mkl)
endif

endif

end subroutine conv_coo2csr
!==========================================================================================
subroutine conv_coo2csr_dp(mat_coo, mat_csr, method)

implicit none

     type (triplet_dp) :: mat_coo
     type (csr_dp) :: mat_csr
     character :: method

     integer,dimension(:),allocatable :: iwork
     integer :: info_mkl, job(8), istat

     call init_sparse_0(mat_csr, mat_coo%ni, mat_coo%nj, mat_coo%nz)

if (method == skit) then

!sparskit convert COO to CSR:
     call coocsr_dp(mat_coo%ni,mat_coo%nz,mat_coo%vals,mat_coo%indi, mat_coo%indj, &
                       mat_csr%vals,mat_csr%indj,mat_csr%indi)
!sparskit sort CSR with column index incr within each row:
     allocate(iwork(2*mat_csr%nz),stat=istat)
     call csort(mat_csr%ni,mat_csr%vals,mat_csr%indj,mat_csr%indi,iwork,.true.)

else

! MKL convert triplet/COO to CSR format
job(1)=2 ! convert COO to CSR and sort
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo%nz ! # of non-zeros in matrix A
job(6)=0 ! Fill all CSR arrays

     call mkl_dcsrcoo(job, mat_coo%ni, mat_csr%vals, mat_csr%indj, mat_csr%indi, mat_csr%nz, &
                 mat_coo%vals, mat_coo%indi, mat_coo%indj, info_mkl)

endif

end subroutine conv_coo2csr_dp

!==========================================================================================
!==========================================================================================
! strange problem with mkl option
subroutine conv_coo2csr_cmplx(mat_coo, mat_csr_cmplx, method)

implicit none

     type (triplet) :: mat_coo
     type (csr_cmplx) :: mat_csr_cmplx
     character :: method

     integer,dimension(:),allocatable :: iwork
     integer :: info_mkl, job(8), istat

     call init_sparse_0(mat_csr_cmplx, mat_coo%ni, mat_coo%nj, mat_coo%nz)

if (method == skit) then

!sparskit convert COO to CSR:
     call coocsr_cmplx(mat_coo%ni,mat_coo%nz,mat_coo%vals,mat_coo%indi, mat_coo%indj, &
                       mat_csr_cmplx%vals,mat_csr_cmplx%indj,mat_csr_cmplx%indi)
!sparskit sort CSR with column index incr within each row:
     allocate(iwork(2*mat_csr_cmplx%nz),stat=istat)
     call csort_cmplx(mat_csr_cmplx%ni,mat_csr_cmplx%vals,mat_csr_cmplx%indj,mat_csr_cmplx%indi,iwork,.true.)

else

! MKL convert triplet/COO to CSR format
job(1)=2 ! convert COO to CSR and sort
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo%nz ! # of non-zeros in matrix A
job(6)=0 ! Fill all CSR arrays

if (wp==dp) then
     call mkl_zcsrcoo(job, mat_coo%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, mat_csr_cmplx%nz, &
                 cmplx(mat_coo%vals, 0., cwp), mat_coo%indi, mat_coo%indj, info_mkl)
elseif (cwp==sp) then
     call mkl_ccsrcoo(job, mat_coo%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, mat_csr_cmplx%nz, &
                 cmplx(mat_coo%vals, 0., cwp), mat_coo%indi, mat_coo%indj, info_mkl)
endif

endif

end subroutine conv_coo2csr_cmplx
!==========================================================================================
subroutine conv_coo2csr_cmplx_cmplx(mat_coo_cmplx, mat_csr_cmplx, method)

implicit none

     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr_cmplx) :: mat_csr_cmplx
     character :: method

     integer,dimension(:),allocatable :: iwork
     integer :: info_mkl, job(8), istat

     call init_sparse_0(mat_csr_cmplx, mat_coo_cmplx%ni, mat_coo_cmplx%nj, mat_coo_cmplx%nz)

if (method == skit) then

!sparskit convert COO to CSR:
     call coocsr_cmplx(mat_coo_cmplx%ni,mat_coo_cmplx%nz,mat_coo_cmplx%vals,mat_coo_cmplx%indi, mat_coo_cmplx%indj, &
                       mat_csr_cmplx%vals,mat_csr_cmplx%indj,mat_csr_cmplx%indi)
!sparskit sort CSR with column index incr within each row:
     allocate(iwork(2*mat_csr_cmplx%nz),stat=istat)
     call csort_cmplx(mat_csr_cmplx%ni,mat_csr_cmplx%vals,mat_csr_cmplx%indj,mat_csr_cmplx%indi,iwork,.true.)

else

! MKL convert triplet/COO to CSR format
job(1)=2 ! convert COO to CSR and sort
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo_cmplx%nz ! # of non-zeros in matrix A
job(6)=0 ! Fill all CSR arrays

if (wp==dp) then
     call mkl_zcsrcoo(job, mat_coo_cmplx%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, mat_csr_cmplx%nz, &
                 mat_coo_cmplx%vals, mat_coo_cmplx%indi, mat_coo_cmplx%indj, info_mkl)
elseif (cwp==sp) then
     call mkl_ccsrcoo(job, mat_coo_cmplx%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, mat_csr_cmplx%nz, &
                 mat_coo_cmplx%vals, mat_coo_cmplx%indi, mat_coo_cmplx%indj, info_mkl)
endif

endif

end subroutine conv_coo2csr_cmplx_cmplx
!==========================================================================================
!==========================================================================================
subroutine conv_csr2coo(mat_csr, mat_coo, method)

implicit none

     type (triplet) :: mat_coo
     type (csr) :: mat_csr
     character :: method

     integer :: info_mkl, job(8), ierr

     call init_sparse_0(mat_coo, mat_csr%ni, mat_csr%nj, mat_csr%nz)

if (method == skit) then
! sparskit convert CSR to COO format
     call csrcoo(mat_csr%ni, 3, mat_csr%nz,mat_csr%vals,mat_csr%indj, mat_csr%indi, &
                       mat_coo%nz, mat_coo%vals,mat_coo%indi,mat_coo%indj, ierr)
else
! MKL convert CSR to triplet/COO format
job(1)=0 ! convert CSR to COO
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo%nz ! # of non-zeros in matrix A
job(6)=3 ! all arrays rowind, colind, acoo are filled in for the output storage.

if (wp==dp) then
     call mkl_dcsrcoo(job, mat_csr%ni, mat_csr%vals, mat_csr%indj, mat_csr%indi, mat_csr%nz, &
                      mat_coo%vals, mat_coo%indi, mat_coo%indj, info_mkl)
elseif (wp==sp) then
     call mkl_scsrcoo(job, mat_csr%ni, mat_csr%vals, mat_csr%indj, mat_csr%indi, mat_csr%nz, &
                      mat_coo%vals, mat_coo%indi, mat_coo%indj, info_mkl)
endif

endif

end subroutine conv_csr2coo
!==========================================================================================
!==========================================================================================
subroutine conv_csr2coo_dp(mat_csr, mat_coo)

implicit none

     type (triplet_dp) :: mat_coo
     type (csr_dp) :: mat_csr
     character :: method

     integer :: info_mkl, job(8), ierr

     call init_sparse_0(mat_coo, mat_csr%ni, mat_csr%nj, mat_csr%nz)

!if (method == skit) then
!! sparskit convert CSR to COO format
!     call csrcoo(mat_csr%ni, 3, mat_csr%nz,mat_csr%vals,mat_csr%indj, mat_csr%indi, &
!                       mat_coo%nz, mat_coo%vals,mat_coo%indi,mat_coo%indj, ierr)
!else
! MKL convert CSR to triplet/COO format
job(1)=0 ! convert CSR to COO
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo%nz ! # of non-zeros in matrix A
job(6)=3 ! all arrays rowind, colind, acoo are filled in for the output storage.

     call mkl_dcsrcoo(job, mat_csr%ni, mat_csr%vals, mat_csr%indj, mat_csr%indi, mat_csr%nz, &
                      mat_coo%vals, mat_coo%indi, mat_coo%indj, info_mkl)
!endif

end subroutine conv_csr2coo_dp
!==========================================================================================
!==========================================================================================
subroutine conv_csr2coo_cmplx(mat_csr_cmplx, mat_coo_cmplx, method)

implicit none

     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr_cmplx) :: mat_csr_cmplx
     character :: method

     integer :: info_mkl, job(8), ierr

     call init_sparse_0(mat_coo_cmplx, mat_csr_cmplx%ni, mat_csr_cmplx%nj, mat_csr_cmplx%nz)

if (method == skit) then
! sparskit convert CSR to COO format
     call csrcoo_cmplx(mat_csr_cmplx%ni, 3, mat_csr_cmplx%nz,mat_csr_cmplx%vals,mat_csr_cmplx%indj, mat_csr_cmplx%indi, &
                       mat_coo_cmplx%nz, mat_coo_cmplx%vals,mat_coo_cmplx%indi,mat_coo_cmplx%indj, ierr)
else
! MKL convert CSR to triplet/COO format
job(1)=0 ! convert CSR to COO
job(2)=1 ! 1-based indexing in CSR
job(3)=1 ! 1-based indexing in COO
job(5)=mat_coo_cmplx%nz ! # of non-zeros in matrix A
job(6)=3 ! all arrays rowind, colind, acoo are filled in for the output storage.

if (cwp==dp) then
     call mkl_zcsrcoo(job, mat_csr_cmplx%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, mat_csr_cmplx%nz, &
                      mat_coo_cmplx%vals, mat_coo_cmplx%indi, mat_coo_cmplx%indj, info_mkl)
elseif (cwp==sp) then
     call mkl_ccsrcoo(job, mat_csr_cmplx%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, mat_csr_cmplx%nz, &
                      mat_coo_cmplx%vals, mat_coo_cmplx%indi, mat_coo_cmplx%indj, info_mkl)
endif

endif

end subroutine conv_csr2coo_cmplx
!==========================================================================================
!==========================================================================================
subroutine conv_csr2csc_cmplx(mat_csr_cmplx, mat_csc_cmplx, method)

! IMPORTANT: mat_csc_cmplx is converted to 0-based for umfpack

implicit none

     type (csc_cmplx) :: mat_csc_cmplx
     type (csr_cmplx) :: mat_csr_cmplx
     character :: method

     integer :: info_mkl, job(8)

     call init_sparse_0(mat_csc_cmplx, mat_csr_cmplx%ni, mat_csr_cmplx%nj, mat_csr_cmplx%nz)

if (method == skit) then
! sparskit convert CSR to CSC format
      call csrcsc_cmplx (mat_csr_cmplx%ni,1,1,mat_csr_cmplx%vals,mat_csr_cmplx%indj,mat_csr_cmplx%indi,&
                   mat_csc_cmplx%vals,mat_csc_cmplx%indi,mat_csc_cmplx%indj)
! convert from 1-based to 0-based
      mat_csc_cmplx%indi = mat_csc_cmplx%indi - 1
      mat_csc_cmplx%indj = mat_csc_cmplx%indj - 1
else
! MKL convert CSR to CSC format
job(1)=0 ! convert CSR to CSC
job(2)=1 ! 1-based indexing in CSR
job(3)=0 ! 1-based indexing in CSC
!job(5)=mat_csc_cmplx%nz ! # of non-zeros in matrix A
job(6)=1 ! all output arrays are filled in for the output storage.

if (wp==dp) then
     call mkl_zcsrcsc(job, mat_csr_cmplx%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, &
                      mat_csc_cmplx%vals,  mat_csc_cmplx%indi, mat_csc_cmplx%indj, info_mkl)
elseif (cwp==sp) then
     call mkl_ccsrcsc(job, mat_csr_cmplx%ni, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, &
                      mat_csc_cmplx%vals,  mat_csc_cmplx%indi, mat_csc_cmplx%indj, info_mkl)
endif

endif

end subroutine conv_csr2csc_cmplx
!==========================================================================================

!==========================================================================================
!------------------- Transform complex sparse NxN into Real sparse 2Nx2N--------------------
!==========================================================================================
subroutine conv_csr_cmplx2csr_dp(csr, a_re)

        type(csr_dp), intent(out)		:: a_re ! a_re is a (2n x 2n) matrix, a_re = [Im(csr) Re(csr); Re(csr) -Im(csr)]
        type (csr_cmplx), intent(in)    :: csr

        type(csr_dp)		       :: csr_re, csr_im
        type(triplet_dp)		   :: coo_re, coo_im, a_re_coo

        integer :: filled, nnz ! points how many elements out of nz are already filled

		call init_sparse_vals(csr_re,csr%indi,csr%indj,real(csr%vals, dp),csr%ni,csr%nj,csr%nz)
		call init_sparse_vals(csr_im,csr%indi,csr%indj,real(imag(csr%vals),dp),csr%ni,csr%nj,csr%nz)
		call dealloc_sparse(csr)

		call csr2coo(csr_re, coo_re)
		call csr2coo(csr_im, coo_im)
		call dealloc_sparse(csr_re)
		call dealloc_sparse(csr_im)

!		Is necessary to remove zero elements in csr_re and csr_im
		call init_sparse_0( a_re_coo, 2*coo_re%ni, 2*coo_re%nj, 2*(count(coo_re%vals .ne. 0.) + count(coo_im%vals .ne. 0.)) )
		! a_re = [Re(csr) -Im(csr); Im(csr) Re(csr)]
		filled=0
		nnz=0
		call concatenate_sparse_dp_rm_nonzero(a_re_coo, coo_im, 0, 0, filled, nnz)
		a_re_coo%indi(filled+1:filled+nnz) = a_re_coo%indi(filled+1-nnz:filled) + coo_re%ni
		a_re_coo%indj(filled+1:filled+nnz) = a_re_coo%indj(filled+1-nnz:filled) + coo_re%nj
		a_re_coo%vals(filled+1:filled+nnz) = -a_re_coo%vals(filled+1-nnz:filled)
		filled = filled + nnz

		call concatenate_sparse_dp_rm_nonzero(a_re_coo, coo_re, 0, coo_re%nj, filled, nnz)
		a_re_coo%indi(filled+1:filled+nnz) = a_re_coo%indi(filled+1-nnz:filled) + coo_re%ni
		a_re_coo%indj(filled+1:filled+nnz) = a_re_coo%indj(filled+1-nnz:filled) - coo_re%nj
		a_re_coo%vals(filled+1:filled+nnz) = a_re_coo%vals(filled+1-nnz:filled)
		filled = filled + nnz

		call dealloc_sparse(coo_im)
		call dealloc_sparse(coo_re)

	!     Save a_re_coo
!	      call save_sparse(a_re_coo, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/'//&
!	      							  'baro_fd/out/0000_00_00__00_00/global/mats/gmres.dat')
!	      print *, "Saved a_re_coo with nnz", a_re_coo%nz

    ! Convert mat from COO to CSR:
	    call coo2csr(a_re_coo, a_re, mkl) 
		call dealloc_sparse(a_re_coo)


end subroutine conv_csr_cmplx2csr_dp
!==========================================================================================

!==========================================================================================
subroutine mul_dia_csr(diag, a_csr, method, b_csr)
! B = Diag * A
implicit none

     type (csr)           :: a_csr, tmp_csr
     type (csr), optional :: b_csr
     real(wp)             :: diag(:)
     character :: method

     integer :: request, sort, info_mkl, j
     character  trans
     integer, allocatable :: tmp(:)

if (method == skit) then

!sparskit: mat=bcmat*mat
    if (present(b_csr)) then
        call init_sparse_0(b_csr, a_csr%ni, a_csr%nj, a_csr%nz)
        call diamua(a_csr%ni, 1, a_csr%vals, a_csr%indj, a_csr%indi, real(diag, kind=wp), &
                                       b_csr%vals, b_csr%indj, b_csr%indi)
    else
       call diamua(a_csr%ni, 0, a_csr%vals, a_csr%indj, a_csr%indi, real(diag, kind=wp), &
                                       a_csr%vals, a_csr%indj, a_csr%indi)
    endif
else

!mkl: mat=bcmat*mat
trans = 'n'
request = 0
sort = 0

     allocate(tmp(a_csr%nj))
     tmp = (/ ( (j), j=1,a_csr%nj) /)

if (wp==dp) then
    if (present(b_csr)) then
          call init_sparse_0(b_csr, a_csr%ni, a_csr%nj, a_csr%nz)
          call mkl_dcsrmultcsr(trans, request, sort, a_csr%ni, a_csr%nj, a_csr%nj, &
                     diag, tmp, (/ tmp, a_csr%nj+1/), &
                     a_csr%vals, a_csr%indj, a_csr%indi, &
                     b_csr%vals, b_csr%indj, b_csr%indi, a_csr%nz, info_mkl)
    else
          call init_sparse_0(tmp_csr, a_csr%ni, a_csr%nj, a_csr%nz)
          call mkl_dcsrmultcsr(trans, request, sort, a_csr%ni, a_csr%nj, a_csr%nj, &
                     diag, tmp, (/ tmp, a_csr%nj+1/), &
                     a_csr%vals, a_csr%indj, a_csr%indi, &
                     tmp_csr%vals, tmp_csr%indj, tmp_csr%indi, tmp_csr%nz, info_mkl)
! copy result back to a_csr
          a_csr = tmp_csr
    endif
elseif (wp==sp) then
    if (present(b_csr)) then
          call init_sparse_0(b_csr, a_csr%ni, a_csr%nj, a_csr%nz)
          call mkl_scsrmultcsr(trans, request, sort, a_csr%ni, a_csr%nj, a_csr%nj, &
                     diag, tmp, (/ tmp, a_csr%nj+1/), &
                     a_csr%vals, a_csr%indj, a_csr%indi, &
                     b_csr%vals, b_csr%indj, b_csr%indi, a_csr%nz, info_mkl)
    else
          call init_sparse_0(tmp_csr, a_csr%ni, a_csr%nj, a_csr%nz)
          call mkl_scsrmultcsr(trans, request, sort, a_csr%ni, a_csr%nj, a_csr%nj, &
                     diag, tmp, (/ tmp, a_csr%nj+1/), &
                     a_csr%vals, a_csr%indj, a_csr%indi, &
                     tmp_csr%vals, tmp_csr%indj, tmp_csr%indi, tmp_csr%nz, info_mkl)
! copy result back to a_csr
          a_csr = tmp_csr
    endif
endif

endif

end subroutine mul_dia_csr
!==========================================================================================
!==========================================================================================
subroutine mul_dia_csr_cmplx(diag, a_csr_cmplx, method, b_csr_cmplx)
! B = Diag * A
implicit none

     type (csr_cmplx)           :: a_csr_cmplx, tmp_csr_cmplx
     type (csr_cmplx), optional :: b_csr_cmplx
     real(wp)                   :: diag(:)
     character :: method

     integer :: request, sort, info_mkl, j
     character  trans
     integer, allocatable :: ones(:)

if (method == skit) then

!sparskit: mat=bcmat*mat
    if (present(b_csr_cmplx)) then
        call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
        call diamua_cmplx(a_csr_cmplx%ni, 1, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, real(diag, kind=wp), &
                                       b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi)
    else
       call diamua_cmplx(a_csr_cmplx%ni, 0, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, real(diag, kind=wp), &
                                       a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi)
    endif
else

!mkl: mat=bcmat*mat
trans = 'n'
request = 0
sort = 0

     allocate(ones(a_csr_cmplx%nj))
     ones = (/ ( (j), j=1,a_csr_cmplx%nj) /)

if (cwp==dp) then
    if (present(b_csr_cmplx)) then
          call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_zcsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, a_csr_cmplx%nz, info_mkl)
    else
          call init_sparse_0(tmp_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_zcsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, tmp_csr_cmplx%indi, tmp_csr_cmplx%nz, info_mkl)
! copy result back to a_csr_cmplx
          a_csr_cmplx = tmp_csr_cmplx
    endif
elseif (cwp==sp) then
    if (present(b_csr_cmplx)) then
          call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_ccsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, a_csr_cmplx%nz, info_mkl)
    else
          call init_sparse_0(tmp_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_ccsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, tmp_csr_cmplx%indi, tmp_csr_cmplx%nz, info_mkl)
! copy result back to a_csr_cmplx
          a_csr_cmplx = tmp_csr_cmplx
    endif
endif

endif

end subroutine mul_dia_csr_cmplx
!==========================================================================================
!==========================================================================================
subroutine mul_dia_int_csr_cmplx(diag, a_csr_cmplx, method, b_csr_cmplx)
! B = Diag * A
implicit none

     type (csr_cmplx)           :: a_csr_cmplx, tmp_csr_cmplx
     type (csr_cmplx), optional :: b_csr_cmplx
     integer                   :: diag(:)
     character :: method

     integer :: request, sort, info_mkl, j
     character  trans
     integer, allocatable :: ones(:)

if (method == skit) then

!sparskit: mat=bcmat*mat
    if (present(b_csr_cmplx)) then
        call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
        call dia_int_mua_cmplx(a_csr_cmplx%ni, 1, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, diag, &
                                       b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi)
    else
       call dia_int_mua_cmplx(a_csr_cmplx%ni, 0, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, diag, &
                                       a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi)
    endif
else

!mkl: mat=bcmat*mat
trans = 'n'
request = 0
sort = 0

     allocate(ones(a_csr_cmplx%nj))
     ones = (/ ( (j), j=1,a_csr_cmplx%nj) /)

if (cwp==dp) then
    if (present(b_csr_cmplx)) then
          call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_zcsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, a_csr_cmplx%nz, info_mkl)
    else
          call init_sparse_0(tmp_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_zcsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, tmp_csr_cmplx%indi, tmp_csr_cmplx%nz, info_mkl)
! copy result back to a_csr_cmplx
          a_csr_cmplx = tmp_csr_cmplx
    endif
elseif (cwp==sp) then
    if (present(b_csr_cmplx)) then
          call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_ccsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, a_csr_cmplx%nz, info_mkl)
    else
          call init_sparse_0(tmp_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
          call mkl_ccsrmultcsr(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nj, &
                     cmplx(diag, 0., cwp), ones, (/ ones, a_csr_cmplx%nj+1/), &
                     a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                     tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, tmp_csr_cmplx%indi, tmp_csr_cmplx%nz, info_mkl)
! copy result back to a_csr_cmplx
          a_csr_cmplx = tmp_csr_cmplx
    endif
endif

endif

end subroutine mul_dia_int_csr_cmplx
!==========================================================================================
!==========================================================================================
subroutine add_csr(a_csr, beta, b_csr, method, c_csr)
! C := A + beta*B
implicit none

     type (csr)           :: a_csr, b_csr, tmp_csr
     type (csr), optional :: c_csr
     real(wp)             :: beta
     character :: method

     integer   :: request, sort, info_mkl, istat
     character :: trans
     integer, allocatable :: iw(:)

if (method == skit) then

!sparskit
    if (present(c_csr)) then
! requires careful pre-allocation, tedious, not done properly
!     call aplsca (mat_csr%ni, mat_csr%nz, mat_csr%vals, mat_csr%indj, mat_csr%indi, &
!                                      - i*pars%omega0, iw)
    else

    endif
else

!mkl
trans = 'n'
request = 1 !the routine computes only values of the array indi of length ni + 1,
            !the memory for this array must be allocated beforehand. On exit the value indi(ni+1) - 1 is nz
sort = 0 ! no reordering

!print "allocated?", allocated(tmp_csr%indi), allocated(tmp_csr%indj), allocated(tmp_csr%vals)


allocate(iw(a_csr%ni + 1),stat=istat)
iw = 0
call init_sparse_0(tmp_csr, 1, 1, 1) ! stupid unnecessary preallocation but sometimes the mkl routine misbehaves

if (wp==dp) then
	     call mkl_dcsradd(trans, request, sort, a_csr%ni, a_csr%nj, &
	                  a_csr%vals, a_csr%indj, a_csr%indi, &
	                  beta, b_csr%vals, b_csr%indj, b_csr%indi, &
	                  tmp_csr%vals, tmp_csr%indj, iw, tmp_csr%nz, info_mkl)

	     call init_sparse_0( tmp_csr, a_csr%ni, a_csr%nj, iw(a_csr%ni+1)-1 )
	     tmp_csr%indi = iw

	deallocate(iw)

	request = 2 !been called previously with the parameter request=1, the output arrays are allocated
	     call mkl_dcsradd(trans, request, sort, a_csr%ni, a_csr%nj, &
	                  a_csr%vals, a_csr%indj, a_csr%indi, &
	                  beta, b_csr%vals, b_csr%indj, b_csr%indi, &
	                  tmp_csr%vals, tmp_csr%indj, tmp_csr%indi, tmp_csr%nz, info_mkl)
elseif (wp==sp) then
	     call mkl_scsradd(trans, request, sort, a_csr%ni, a_csr%nj, &
	                  a_csr%vals, a_csr%indj, a_csr%indi, &
	                  beta, b_csr%vals, b_csr%indj, b_csr%indi, &
	                  tmp_csr%vals, tmp_csr%indj, iw, tmp_csr%nz, info_mkl)

	     call init_sparse_0( tmp_csr, a_csr%ni, a_csr%nj, iw(a_csr%ni+1)-1 )
	     tmp_csr%indi = iw

	deallocate(iw)

	request = 2 !been called previously with the parameter request=1, the output arrays are allocated
	     call mkl_scsradd(trans, request, sort, a_csr%ni, a_csr%nj, &
	                  a_csr%vals, a_csr%indj, a_csr%indi, &
	                  beta, b_csr%vals, b_csr%indj, b_csr%indi, &
	                  tmp_csr%vals, tmp_csr%indj, tmp_csr%indi, tmp_csr%nz, info_mkl)
endif

    if (present(c_csr)) then
          call init_sparse_0(c_csr, tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
! copy result to c_csr
          c_csr = tmp_csr
    else
          call dealloc_sparse(a_csr)
          call init_sparse_0(a_csr, tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
! copy result back to a_csr
          a_csr = tmp_csr
    endif

endif

end subroutine add_csr
!==========================================================================================
!==========================================================================================
subroutine add_csr_cmplx(a_csr_cmplx, beta, b_csr_cmplx, method, c_csr_cmplx)
! C := A + beta*B
implicit none

     type (csr_cmplx)           :: a_csr_cmplx, b_csr_cmplx, tmp_csr_cmplx
     type (csr_cmplx), optional :: c_csr_cmplx
     complex(cwp)               :: beta
     character :: method

     integer   :: request, sort, info_mkl, istat
     character :: trans
     integer, allocatable :: iw(:)

if (method == skit) then

!sparskit
    if (present(c_csr_cmplx)) then
! requires careful pre-allocation, tedious, not done properly
!     call aplsca_cmplx (mat_csr_cmplx%ni, mat_csr_cmplx%nz, mat_csr_cmplx%vals, mat_csr_cmplx%indj, mat_csr_cmplx%indi, &
!                                      - i*pars%omega0, iw)
    else

    endif
else

!mkl
trans = 'n'
request = 1 !the routine computes only values of the array indi of length ni + 1,
            !the memory for this array must be allocated beforehand. On exit the value indi(ni+1) - 1 is nz
sort = 0 ! no reordering

!print "allocated?", allocated(tmp_csr_cmplx%indi), allocated(tmp_csr_cmplx%indj), allocated(tmp_csr_cmplx%vals)


allocate(iw(a_csr_cmplx%ni + 1),stat=istat)
iw = 0
call init_sparse_0(tmp_csr_cmplx, 1, 1, 1) ! stupid unnecessary preallocation but sometimes the mkl routine misbehaves

if (cwp==dp) then
	     call mkl_zcsradd(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, &
	                  a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
	                  beta, b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, &
	                  tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, iw, tmp_csr_cmplx%nz, info_mkl)

	     call init_sparse_0( tmp_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, iw(a_csr_cmplx%ni+1)-1 )
	     tmp_csr_cmplx%indi = iw

	deallocate(iw)

	request = 2 !been called previously with the parameter request=1, the output arrays are allocated
	     call mkl_zcsradd(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, &
	                  a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
	                  beta, b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, &
	                  tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, tmp_csr_cmplx%indi, tmp_csr_cmplx%nz, info_mkl)
elseif (cwp==sp) then
	     call mkl_ccsradd(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, &
	                  a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
	                  beta, b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, &
	                  tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, iw, tmp_csr_cmplx%nz, info_mkl)

	     call init_sparse_0( tmp_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, iw(a_csr_cmplx%ni+1)-1 )
	     tmp_csr_cmplx%indi = iw

	deallocate(iw)

	request = 2 !been called previously with the parameter request=1, the output arrays are allocated
	     call mkl_ccsradd(trans, request, sort, a_csr_cmplx%ni, a_csr_cmplx%nj, &
	                  a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
	                  beta, b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, &
	                  tmp_csr_cmplx%vals, tmp_csr_cmplx%indj, tmp_csr_cmplx%indi, tmp_csr_cmplx%nz, info_mkl)
endif

    if (present(c_csr_cmplx)) then
          call init_sparse_0(c_csr_cmplx, tmp_csr_cmplx%ni, tmp_csr_cmplx%nj, tmp_csr_cmplx%nz)
! copy result to c_csr_cmplx
          c_csr_cmplx = tmp_csr_cmplx
    else
          call dealloc_sparse(a_csr_cmplx)
          call init_sparse_0(a_csr_cmplx, tmp_csr_cmplx%ni, tmp_csr_cmplx%nj, tmp_csr_cmplx%nz)
! copy result back to a_csr_cmplx
          a_csr_cmplx = tmp_csr_cmplx
    endif

endif

end subroutine add_csr_cmplx
!==========================================================================================
!==========================================================================================
subroutine add_diacsr_cmplx(a_csr_cmplx, diag, b_csr_cmplx)
! NO MEMORY ALLOCATION, diagonal elements in A are nonzeros
! B = A + Diag
!or
! A = A + Diag
implicit none

     type (csr_cmplx)              :: a_csr_cmplx
     type (csr_cmplx), optional    :: b_csr_cmplx
     complex(cwp)                  :: diag(:)
!     character :: method

!     integer :: request, sort, info_mkl, j
!     character  trans
!     integer, allocatable :: ones(:)

     integer :: iw(a_csr_cmplx%ni)


!if (method == skit) then

!sparskit
    if (present(b_csr_cmplx)) then
        call init_sparse_0(b_csr_cmplx, a_csr_cmplx%ni, a_csr_cmplx%nj, a_csr_cmplx%nz)
!        print *, "2"
        call apldia_cmplx(a_csr_cmplx%ni, 1, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                          diag, b_csr_cmplx%vals, b_csr_cmplx%indj, b_csr_cmplx%indi, iw)
    else
!    print *, "1"
        call apldia_cmplx(a_csr_cmplx%ni, 0, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, &
                          diag, a_csr_cmplx%vals, a_csr_cmplx%indj, a_csr_cmplx%indi, iw)
    endif
!else
!
!mkl
!
!endif

end subroutine add_diacsr_cmplx
!==========================================================================================
!==========================================================================================

end module my_sparse_aggregate
