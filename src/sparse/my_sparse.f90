module my_sparse

     use precisions, only: sp, dp, wp, cwp
	 use dispmodule

     character, parameter :: mkl = 'm' ! MKL functions and routines
     character, parameter :: skit = 's' ! SparseKit and ad-hoc functions

!==========================================================================================
     type triplet_int
!          integer(kind=LongIntType)  :: ni, nj, nz
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          real(wp),allocatable         :: vals(:)
     end type

     type triplet
!          integer(kind=LongIntType)  :: ni, nj, nz
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          real(wp),allocatable         :: vals(:)
     end type

     type triplet_sp
!          integer(kind=LongIntType)  :: ni, nj, nz
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          real(sp),allocatable         :: vals(:)
     end type

     type triplet_dp
!          integer(kind=LongIntType)  :: ni, nj, nz
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          real(dp),allocatable         :: vals(:)
     end type

     type triplet_cmplx
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          complex(cwp),allocatable      :: vals(:)
     end type

     type csr
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          real(wp),allocatable      :: vals(:)
     end type

     type csr_dp
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          real(dp),allocatable      :: vals(:)
     end type

     type csr_cmplx
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          complex(cwp),allocatable      :: vals(:)
     end type

     type csc_cmplx
          integer                    :: ni, nj, nz
          integer,allocatable        :: indi(:), indj(:)
          complex(cwp),allocatable   :: vals(:)
     end type
!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!---------------------------- INITIATE SPARSE ---------------------------------------
     interface init_sparse_0
          module procedure csr_init_0
          module procedure csr_dp_init_0
          module procedure csr_cmplx_init_0
          module procedure csc_cmplx_init_0
          module procedure init_sparse_zeros
          module procedure init_sparse_zeros_sp
          module procedure init_sparse_zeros_dp
          module procedure init_sparse_zeros_int
          module procedure init_sparse_zeros_cmplx
     end interface  
!---------------------------- INITIATE WITH VALS -------------------------------------
     interface init_sparse_vals
          module procedure init_sparse
          module procedure init_sparse_cmplx
          module procedure init_sparse_csr_dp
          module procedure init_sparse_csr_cmplx
          module procedure init_sparse_csr_cmplx_cmplx
     end interface
!---------------------------- INITIATE DIAGONAL -------------------------------------
     interface init_speye
          module procedure speye_csr_cmplx
          module procedure speye_csr_vec_cmplx
     end interface
!---------------------------- DEALLOCATE SPARSE -------------------------------------
     interface dealloc_sparse
          module procedure deallocate_sparse
          module procedure deallocate_sparse_dp
          module procedure deallocate_sparse_int
          module procedure deallocate_sparse_cmplx
          module procedure deallocate_csr_cmplx
          module procedure deallocate_csr
          module procedure deallocate_csr_dp
          module procedure deallocate_csc_cmplx
     end interface
!---------------------------- TRIM SPARSE -------------------------------------
!     interface trim_sparse
!          module procedure trim_sparse_dp
!     end interface
!---------------------------- TRIPLET VEC MULT -------------------------------
     interface coo_vec_mul
          module procedure mul_coo_vec
          module procedure mul_coo_sp_vec_sp
          module procedure mul_coo_dp_vec_dp
          module procedure mul_coo_vec_cmplx
          module procedure mul_coo_vec_cmplx_cmplx
      end interface
!---------------------------- VEC TRIPLET MULT -------------------------------
     interface vec_coo_mul
          module procedure vector_triplet_int_mult
          module procedure vector_cmplx_triplet_int_mult
      end interface
!---------------------------- TRIPLET MATRIX MULT -----------------------------------
     interface coo_mat_mul
		module procedure mul_coo_int_mat
		module procedure mul_coo_int_mat_cmplx
      end interface
!---------------------------- MATRIX TRIPLET MULT -----------------------------------
     interface mat_coo_mul
		module procedure mul_mat_coo_int
		module procedure mul_mat_cmplx_coo_int
      end interface
!==========================================================================================

contains

!==========================================================================================

subroutine csr_init_0(A, ni,nj, nz)

          integer :: istat, j
          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (csr), intent(inout)            :: A ! matrix in csr form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(ni+1), stat = istat)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = istat)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = istat)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi = nz + 1
A%indi(1:(nz-1)/nj+1) = (/ ( ((j-1)*nj + 1), j=1,(nz-1)/nj+1 ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine csr_init_0
!==========================================================================================
!==========================================================================================

subroutine csr_dp_init_0(A, ni,nj, nz)

          integer :: istat, j
          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (csr_dp), intent(inout)            :: A ! matrix in csr form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(ni+1), stat = istat)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = istat)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = istat)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi = nz + 1
A%indi(1:(nz-1)/nj+1) = (/ ( ((j-1)*nj + 1), j=1,(nz-1)/nj+1 ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine csr_dp_init_0
!==========================================================================================
!==========================================================================================

subroutine csr_cmplx_init_0(A, ni,nj, nz)

          integer :: istat, j
          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (csr_cmplx), intent(inout)            :: A ! matrix in csr form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(ni+1), stat = istat)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = istat)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = istat)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi = nz + 1
A%indi(1:(nz-1)/nj+1) = (/ ( ((j-1)*nj + 1), j=1,(nz-1)/nj+1 ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine csr_cmplx_init_0
!==========================================================================================
!==========================================================================================

subroutine csc_cmplx_init_0(A, ni,nj, nz)

          integer :: istat, j
          integer                         :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (csc_cmplx), intent(inout) :: A ! matrix in csc form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indj(nj+1), stat = istat)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indi(nz), stat = istat)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = istat)

A%ni = ni
A%nj = nj
A%nz = nz

A%indj = nz + 1
A%indj(1:(nz-1)/ni+1) = (/ ( ((j-1)*ni + 1), j=1,(nz-1)/ni+1 ) /)
A%indi(:) = (/ ( (mod(j-1,ni) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine csc_cmplx_init_0
!==========================================================================================
!==========================================================================================

subroutine init_sparse_zeros(A, ni,nj, nz)

          integer :: statusi, statusj, statusval, j

          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (triplet), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = (/ ( ((j-1)/nj + 1), j=1,nz ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine init_sparse_zeros

subroutine init_sparse_zeros_sp(A, ni,nj, nz)

          integer :: statusi, statusj, statusval, j

          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (triplet_sp), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = (/ ( ((j-1)/nj + 1), j=1,nz ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine init_sparse_zeros_sp


subroutine init_sparse_zeros_dp(A, ni,nj, nz)

          integer :: statusi, statusj, statusval, j

          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (triplet_dp), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = (/ ( ((j-1)/nj + 1), j=1,nz ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine init_sparse_zeros_dp

!==========================================================================================
!==========================================================================================

subroutine init_sparse_zeros_int(A, ni,nj, nz)

          integer :: statusi, statusj, statusval, j

          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (triplet_int), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = (/ ( ((j-1)/nj + 1), j=1,nz ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0

end subroutine init_sparse_zeros_int

!==========================================================================================
!==========================================================================================

subroutine init_sparse_zeros_cmplx(A, ni,nj, nz)

          integer :: statusi, statusj, statusval, j

          integer                  :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          type (triplet_cmplx), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = (/ ( ((j-1)/nj + 1), j=1,nz ) /)
A%indj(:) = (/ ( (mod(j-1,nj) + 1), j=1,nz ) /)
A%vals(:) = 0


end subroutine init_sparse_zeros_cmplx

!==========================================================================================
!==========================================================================================
subroutine deallocate_sparse(A)
          type (triplet) :: A

deallocate(A%indi, A%indj, A%vals)
end subroutine deallocate_sparse

subroutine deallocate_sparse_dp(A)
          type (triplet_dp) :: A

deallocate(A%indi, A%indj, A%vals)
end subroutine deallocate_sparse_dp

subroutine deallocate_sparse_int(A)

          type (triplet_int) :: A

deallocate(A%indi, A%indj, A%vals)

end subroutine deallocate_sparse_int


subroutine deallocate_sparse_cmplx(A)

          type (triplet_cmplx) :: A

deallocate(A%indi, A%indj, A%vals)

end subroutine deallocate_sparse_cmplx

subroutine deallocate_csr(A)

          type (csr) :: A

deallocate(A%indi, A%indj, A%vals)
end subroutine deallocate_csr
subroutine deallocate_csr_dp(A)

          type (csr_dp) :: A

deallocate(A%indi, A%indj, A%vals)
end subroutine deallocate_csr_dp

subroutine deallocate_csr_cmplx(A)

          type (csr_cmplx) :: A

deallocate(A%indi, A%indj, A%vals)

end subroutine deallocate_csr_cmplx

subroutine deallocate_csc_cmplx(A)

          type (csc_cmplx) :: A

deallocate(A%indi, A%indj, A%vals)

end subroutine deallocate_csc_cmplx
!==========================================================================================
!==========================================================================================

subroutine init_sparse(A,indi,indj,vals,ni,nj, nz)

          integer :: statusi, statusj, statusval

          integer          :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          integer          :: indi(:), indj(:) ! i and j index for the elements
          real(wp)         :: vals(:)       ! elements' values
          type (triplet), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = indi(1:nz)
A%indj(:) = indj(1:nz)
A%vals(:) = vals(1:nz)

end subroutine init_sparse

subroutine init_sparse_dp(A,indi,indj,vals,ni,nj, nz)

          integer :: statusi, statusj, statusval

          integer          :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          integer          :: indi(:), indj(:) ! i and j index for the elements
          real(dp)         :: vals(:)       ! elements' values
          type (triplet_dp), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = indi(1:nz)
A%indj(:) = indj(1:nz)
A%vals(:) = vals(1:nz)

end subroutine init_sparse_dp

!==========================================================================================
!==========================================================================================

subroutine init_sparse_cmplx(A,indi,indj,vals,ni,nj, nz)

          integer :: statusi, statusj, statusval

          integer			:: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          integer 			:: indi(:), indj(:) ! i and j index for the elements
          complex(cwp)      :: vals(:)       ! elements' values
          type (triplet_cmplx), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(nz), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = indi(1:nz)
A%indj(:) = indj(1:nz)
A%vals(:) = vals(1:nz)

end subroutine init_sparse_cmplx

!==========================================================================================
!==========================================================================================

subroutine init_sparse_csr_dp(A,indi,indj,vals,ni,nj, nz)

          integer :: statusi, statusj, statusval

          integer             :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          integer             :: indi(:), indj(:) ! i and j index for the elements
          real(dp)            :: vals(:)       ! elements' values
          type (csr_dp), intent(inout)            :: A ! matrix

    allocate(A%indi(ni+1), stat = statusi)
    allocate(A%indj(nz), stat = statusj)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = indi(1:ni+1)
A%indj(:) = indj(1:nz)
A%vals(:) = vals(1:nz)

end subroutine init_sparse_csr_dp

!==========================================================================================
!==========================================================================================

subroutine init_sparse_csr_cmplx(A,indi,indj,vals,ni,nj, nz)

          integer :: statusi, statusj, statusval

          integer             :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          integer             :: indi(:), indj(:) ! i and j index for the elements
          real(wp)            :: vals(:)       ! elements' values
          type (csr_cmplx), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(ni+1), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(nz), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(nz), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = indi(1:ni+1)
A%indj(:) = indj(1:nz)
A%vals(:) = cmplx(vals(1:nz),0., kind=cwp)

end subroutine init_sparse_csr_cmplx

!==========================================================================================
!==========================================================================================

subroutine init_sparse_csr_cmplx_cmplx(A,indi,indj,vals,ni,nj, nz)

          integer :: statusi, statusj, statusval

          integer             :: ni, nj, nz ! ni, nj -- number of rows and columns. nz -- number of elements
          integer             :: indi(:), indj(:) ! i and j index for the elements
          complex(cwp)        :: vals(:)       ! elements' values
          type (csr_cmplx), intent(inout)            :: A ! matrix in triplet form

!    if (associated(A%indi)) deallocate(A%indi)
    allocate(A%indi(size(indi)), stat = statusi)
!    if (associated(A%indj)) deallocate(A%indj)
    allocate(A%indj(size(indj)), stat = statusj)
!    if (associated(A%vals)) deallocate(A%vals)
    allocate(A%vals(size(vals)), stat = statusval)

A%ni = ni
A%nj = nj
A%nz = nz

A%indi(:) = indi(:)
A%indj(:) = indj(:)
A%vals(:) = vals(:)

end subroutine init_sparse_csr_cmplx_cmplx

!==========================================================================================
!==========================================================================================

function speye_csr_cmplx(n,const) ! sparse matrix with const on the diagonal

          integer   :: n, j
          complex(cwp) :: const
          type (csr_cmplx) :: speye_csr_cmplx

call init_sparse_0(speye_csr_cmplx, n,n,n)

speye_csr_cmplx%indi(:) = (/ ( (j), j=1,n+1 ) /)
speye_csr_cmplx%indj(:) = (/ ( (j), j=1,n ) /)
speye_csr_cmplx%vals(:) = const

end function speye_csr_cmplx
!==========================================================================================
function speye_csr_vec_cmplx(n,vec) ! sparse matrix with const on the diagonal

          integer   	:: n, j
          complex(cwp) 	:: vec(n)
          type (csr_cmplx) :: speye_csr_vec_cmplx

call init_sparse_0(speye_csr_vec_cmplx, n,n,n)

speye_csr_vec_cmplx%indi(:) = (/ ( (j), j=1,n+1 ) /)
speye_csr_vec_cmplx%indj(:) = (/ ( (j), j=1,n ) /)
speye_csr_vec_cmplx%vals(:) = vec(:)

end function speye_csr_vec_cmplx

!==========================================================================================
!==========================================================================================
subroutine concatenate_sparse(A, B, shift_i, shift_j, filled)

          integer :: shift_i, shift_j ! shift indices of B when filling out A
          integer :: filled ! points how many elements out of nz are already filled

          type (triplet), intent(inout)            :: A ! bigger matrix
          type (triplet), intent(in)            :: B ! will be a submatrix of A

A%indi(filled+1:filled+B%nz) = B%indi(:) + shift_i
A%indj(filled+1:filled+B%nz) = B%indj(:) + shift_j
A%vals(filled+1:filled+B%nz) = B%vals(:)

filled = filled + B%nz

end subroutine concatenate_sparse
!==========================================================================================
!==========================================================================================
subroutine concatenate_sparse_dp_rm_nonzero(A, B, shift_i, shift_j, filled, nnz)

          integer :: shift_i, shift_j ! shift indices of B when filling out A
          integer :: filled, nnz ! points how many elements out of nz are already filled

          type (triplet_dp), intent(inout)            :: A ! bigger matrix
          type (triplet_dp), intent(in)            :: B ! will be a submatrix of A

          integer :: istat
          integer,allocatable   :: indi(:), indj(:)
          real(dp),allocatable  :: vals(:)

nnz = count(B%vals .ne. 0.)
!print *, nnz, B%nz
    allocate(indi(B%nz),indj(B%nz),vals(B%nz), stat = istat)
    indi = 0
    indj = 0
    vals = 0
where (B%vals .ne. 0)
	indi = B%indi
	indj = B%indj
	vals = B%vals
end where

	A%indi(filled+1:filled+nnz) = pack(indi, indi .ne. 0) + shift_i
	A%indj(filled+1:filled+nnz) = pack(indj, indj .ne. 0) + shift_j
	A%vals(filled+1:filled+nnz) = pack(vals, vals .ne. 0)

filled = filled + nnz

end subroutine concatenate_sparse_dp_rm_nonzero

!==========================================================================================
!==========================================================================================

function scalar_mult(A, b, c)
!    scalar_mult = b*A + c

implicit none

          real(wp), intent(in)          :: b, c

          type (triplet), intent(in)  :: A
          type (triplet)              :: scalar_mult

scalar_mult = A !initialize scalar_mult

scalar_mult%vals(:) = b*A%vals(:) + c

end function scalar_mult

!==========================================================================================

subroutine diapos ( n, nz, ja, ia, idiag )
!*****************************************************************************80
!returns the positions of the diagonal elements of a sparse matrix.

!    Output, integer IDIAG(N); the I-th entry of IDIAG points to the
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element
!    is found, the entry is set to 0.
!
  implicit none

  integer n, nz

  integer j
  integer ia(nz), ja(nz)
  integer idiag(n)

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  do j = 1, nz
      if ( ja(j) == ia(j) ) then
        idiag(ja(j)) = j
      end if
  end do

end subroutine diapos

!==========================================================================================
!==========================================================================================
subroutine right_mult_diag(A, diag, nj)

          integer :: cj
          integer, intent(in)           :: nj
          real(wp), intent(in)            :: diag(nj)

          type (triplet), intent(inout) :: A

if (nj /= A%nj) then
      write(*, '("Error in right_mult_diag", i10, "/=", i10)'), nj,A%nj
      print *,A%ni, A%nj
      stop
endif

  do cj = 1,A%nz
       A%vals(cj) = A%vals(cj) * diag(A%indj(cj))
  end do

end subroutine right_mult_diag

!==========================================================================================
!==========================================================================================

subroutine left_mult_diag(A, diag, ni)

          integer :: ci
          integer, intent(in)           :: ni
          real(wp), intent(in)            :: diag(ni)

          type (triplet), intent(inout) :: A

if (ni /= A%ni) then
      write(*, '("Error in left_mult_diag")')
      stop
endif

  do ci = 1,A%nz
       A%vals(ci) = diag(A%indi(ci)) * A%vals(ci)
  end do

end subroutine left_mult_diag

!==========================================================================================
!==========================================================================================

subroutine left_right_mult_diag(A, diag_l, ni, diag_r, nj)

          integer :: c
          integer, intent(in)           :: ni, nj
          real(wp), intent(in)            :: diag_l(ni), diag_r(nj)

          type (triplet), intent(inout) :: A

if ((ni /= A%ni).or.(nj /= A%nj)) then
      write(*, '("Error in left_right_mult_diag")')
      stop
endif

  do c = 1,A%nz
       A%vals(c) = diag_l(A%indi(c)) * A%vals(c) * diag_r(A%indj(c))
  end do

end subroutine left_right_mult_diag

subroutine left_right_mult_diag_dp(A, diag_l, ni, diag_r, nj)

          integer :: c
          integer, intent(in)           :: ni, nj
          real(dp), intent(in)            :: diag_l(ni), diag_r(nj)

          type (triplet_dp), intent(inout) :: A

if ((ni /= A%ni).or.(nj /= A%nj)) then
      write(*, '("Error in left_right_mult_diag")')
      stop
endif

  do c = 1,A%nz
       A%vals(c) = diag_l(A%indi(c)) * A%vals(c) * diag_r(A%indj(c))
  end do

end subroutine left_right_mult_diag_dp
!==========================================================================================

!==========================================================================================
!---------------------------- TRIPLET-VECTOR MULTIPLY --------------------------
!==========================================================================================
subroutine mul_coo_vec(a_coo, b, method, c)
! C := A*b
!or
! b := A*b
implicit none

     type (triplet)                       :: a_coo
     real(wp)                         :: b(:)
     real(wp), allocatable            :: tmp(:)
     real(wp), optional               :: c(:)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     real(wp):: alpha, beta


if (method == skit) then

!sparskit
    if (present(c)) then
         call triplet_vector_mult(a_coo, b, c)
    else
         allocate(tmp(a_coo%ni), stat = istatus)
         call triplet_vector_mult(a_coo, b, tmp)
          ! copy result back to b
         b = tmp
    endif
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

if (wp==dp) then
    if (present(c)) then
                          !print *, vec_vec_dotproduct(c, c),vec_vec_dotproduct(b, b)
          call mkl_dcoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, c)
    else
          allocate(tmp(a_coo%ni),  stat = istatus)
          call mkl_dcoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, tmp)
          ! copy result back to b
          b = tmp
    endif
elseif (wp == sp) then
    if (present(c)) then
                          !print *, vec_vec_dotproduct(c, c),vec_vec_dotproduct(b, b)
          call mkl_scoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, c)
    else
          allocate(tmp(a_coo%ni),  stat = istatus)
          call mkl_scoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, tmp)
          ! copy result back to b
          b = tmp
    endif
endif

endif

end subroutine mul_coo_vec
!==========================================================================================
    subroutine triplet_vector_mult(A,x,b)
! multiply  A * x = b
  type (triplet), intent(in)        :: A
  real(wp),dimension(A%nj),intent(in)  :: x
  real(wp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_vector_mult
  !==========================================================================================
subroutine mul_coo_sp_vec_sp(a_coo, b, method, c)
! C := A*b
!or
! b := A*b
implicit none

     type (triplet_sp)                       :: a_coo
     real(sp)                         :: b(:)
     real(sp), allocatable            :: tmp(:)
     real(sp), optional               :: c(:)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     real(sp):: alpha, beta


if (method == skit) then

!sparskit
    if (present(c)) then
         call triplet_vector_mult_sp(a_coo, b, c)
    else
         allocate(tmp(a_coo%ni), stat = istatus)
         call triplet_vector_mult_sp(a_coo, b, tmp)
          ! copy result back to b
         b = tmp
    endif
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

    if (present(c)) then
                          !print *, vec_vec_dotproduct(c, c),vec_vec_dotproduct(b, b)
          call mkl_scoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, c)
    else
          allocate(tmp(a_coo%ni),  stat = istatus)
          call mkl_scoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, tmp)
          ! copy result back to b
          b = tmp
    endif

endif

end subroutine mul_coo_sp_vec_sp
!==========================================================================================
subroutine triplet_vector_mult_sp(A,x,b)
! multiply  A * x = b
  type (triplet_sp), intent(in)        :: A
  real(sp),dimension(A%nj),intent(in)  :: x
  real(sp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_vector_mult_sp
!==========================================================================================
subroutine mul_coo_dp_vec_dp(a_coo, b, method, c)
! C := A*b
!or
! b := A*b
implicit none

     type (triplet_dp)                       :: a_coo
     real(dp)                         :: b(:)
     real(dp), allocatable            :: tmp(:)
     real(dp), optional               :: c(:)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     real(dp):: alpha, beta


if (method == skit) then

!sparskit
    if (present(c)) then
         call triplet_vector_mult_dp(a_coo, b, c)
    else
         allocate(tmp(a_coo%ni), stat = istatus)
         call triplet_vector_mult_dp(a_coo, b, tmp)
          ! copy result back to b
         b = tmp
    endif
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

    if (present(c)) then
                          !print *, vec_vec_dotproduct(c, c),vec_vec_dotproduct(b, b)
          call mkl_dcoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, c)
    else
          allocate(tmp(a_coo%ni),  stat = istatus)
          call mkl_dcoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          a_coo%vals, a_coo%indi, a_coo%indj, a_coo%nz, b, beta, tmp)
          ! copy result back to b
          b = tmp
    endif

endif

end subroutine mul_coo_dp_vec_dp
!==========================================================================================
subroutine triplet_vector_mult_dp(A,x,b)
! multiply  A * x = b
  type (triplet_dp), intent(in)        :: A
  real(dp),dimension(A%nj),intent(in)  :: x
  real(dp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_vector_mult_dp

!==========================================================================================
!!!!!!!!!!!!!!!!!!	WARNING	!!!!!!!!!!!!!!!!!!
! 	Same as below but simpler

!subroutine mul_coo_cmplx_vec(a_coo, b, method, c)
!	! C := A*b
!	!or
!	! b := A*b
!implicit none
!
!     type (triplet)                       :: a_coo
!     complex(cwp)                         :: b(:)
!     real(wp), allocatable            :: tmp1(:), tmp2(:)
!     complex(cwp), optional               :: c(:)
!     character :: method
!
!     integer   :: istatus
!
!	allocate(tmp1(a_coo%ni),tmp2(a_coo%ni), stat=istatus)
!	tmp1 = real(b)
!	tmp2 = imag(b)
!
!	call mul_coo_vec(a_coo, tmp1, method)
!	call mul_coo_vec(a_coo, tmp2, method)
!
!	if (present(c)) then
!		c = cmplx(tmp1, tmp2, kind=cwp)
!	else
!		b = cmplx(tmp1, tmp2, kind=cwp)
!	endif
!
!end subroutine mul_coo_cmplx_vec

!==========================================================================================
subroutine mul_coo_vec_cmplx_cmplx(a_coo_cmplx, b_cmplx, method, c_cmplx)
! C := A*b
!or
! b := A*b
implicit none

     type (triplet_cmplx)                       :: a_coo_cmplx
     complex(cwp)                         :: b_cmplx(:)
     complex(cwp), allocatable            :: tmp_cmplx(:)
     complex(cwp), optional               :: c_cmplx(:)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     complex(cwp):: alpha, beta


if (method == skit) then

!sparskit
    if (present(c_cmplx)) then
         call triplet_cmplx_vector_cmplx_mult(a_coo_cmplx, b_cmplx, c_cmplx)
    else
         allocate(tmp_cmplx(a_coo_cmplx%ni), stat = istatus)
         call triplet_cmplx_vector_cmplx_mult(a_coo_cmplx, b_cmplx, tmp_cmplx)
          ! copy result back to b_cmplx
         b_cmplx = tmp_cmplx
    endif
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

if (cwp==dp) then
    if (present(c_cmplx)) then
                          !print *, vec_vec_dotproduct(c_cmplx, c_cmplx),vec_vec_dotproduct(b_cmplx, b_cmplx)
          call mkl_zcoomv(transa, a_coo_cmplx%ni, a_coo_cmplx%nj, alpha, matdescra, &
                          a_coo_cmplx%vals, a_coo_cmplx%indi, a_coo_cmplx%indj, a_coo_cmplx%nz, b_cmplx, beta, c_cmplx)
    else
          allocate(tmp_cmplx(a_coo_cmplx%ni),  stat = istatus)
          call mkl_zcoomv(transa, a_coo_cmplx%ni, a_coo_cmplx%nj, alpha, matdescra, &
                          a_coo_cmplx%vals, a_coo_cmplx%indi, a_coo_cmplx%indj, a_coo_cmplx%nz, b_cmplx, beta, tmp_cmplx)
          ! copy result back to b_cmplx
          b_cmplx = tmp_cmplx
    endif
elseif (cwp == sp) then
    if (present(c_cmplx)) then
                          !print *, vec_vec_dotproduct(c_cmplx, c_cmplx),vec_vec_dotproduct(b_cmplx, b_cmplx)
          call mkl_ccoomv(transa, a_coo_cmplx%ni, a_coo_cmplx%nj, alpha, matdescra, &
                          a_coo_cmplx%vals, a_coo_cmplx%indi, a_coo_cmplx%indj, a_coo_cmplx%nz, b_cmplx, beta, c_cmplx)
    else
          allocate(tmp_cmplx(a_coo_cmplx%ni),  stat = istatus)
          call mkl_ccoomv(transa, a_coo_cmplx%ni, a_coo_cmplx%nj, alpha, matdescra, &
                          a_coo_cmplx%vals, a_coo_cmplx%indi, a_coo_cmplx%indj, a_coo_cmplx%nz, b_cmplx, beta, tmp_cmplx)
          ! copy result back to b_cmplx
          b_cmplx = tmp_cmplx
    endif
endif


endif

end subroutine mul_coo_vec_cmplx_cmplx
!==========================================================================================
subroutine mul_coo_vec_cmplx(a_coo, b_cmplx, method, c_cmplx)
! C := A*b
!or
! b := A*b
implicit none

     type (triplet)                       :: a_coo
     complex(cwp)                         :: b_cmplx(:)
     complex(cwp), allocatable            :: tmp_cmplx(:)
     complex(cwp), optional               :: c_cmplx(:)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     complex(cwp):: alpha, beta


if (method == skit) then

!sparskit
    if (present(c_cmplx)) then
         call triplet_vector_cmplx_mult(a_coo, b_cmplx, c_cmplx)
    else
         allocate(tmp_cmplx(a_coo%ni), stat = istatus)
         call triplet_vector_cmplx_mult(a_coo, b_cmplx, tmp_cmplx)
          ! copy result back to b_cmplx
         b_cmplx = tmp_cmplx
    endif
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

if (cwp==dp) then
    if (present(c_cmplx)) then
                          !print *, vec_vec_dotproduct(c_cmplx, c_cmplx),vec_vec_dotproduct(b_cmplx, b_cmplx)
          call mkl_zcoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          cmplx(a_coo%vals, 0, kind=dp), a_coo%indi, a_coo%indj, a_coo%nz, b_cmplx, beta, c_cmplx)
    else
          allocate(tmp_cmplx(a_coo%ni),  stat = istatus)
          call mkl_zcoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          cmplx(a_coo%vals, 0, kind=dp), a_coo%indi, a_coo%indj, a_coo%nz, b_cmplx, beta, tmp_cmplx)
          ! copy result back to b_cmplx
          b_cmplx = tmp_cmplx
    endif
elseif (cwp == sp) then
    if (present(c_cmplx)) then
                          !print *, vec_vec_dotproduct(c_cmplx, c_cmplx),vec_vec_dotproduct(b_cmplx, b_cmplx)
          call mkl_ccoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          cmplx(a_coo%vals, 0, kind=sp), a_coo%indi, a_coo%indj, a_coo%nz, b_cmplx, beta, c_cmplx)
    else
          allocate(tmp_cmplx(a_coo%ni),  stat = istatus)
          call mkl_ccoomv(transa, a_coo%ni, a_coo%nj, alpha, matdescra, &
                          cmplx(a_coo%vals, 0, kind=sp), a_coo%indi, a_coo%indj, a_coo%nz, b_cmplx, beta, tmp_cmplx)
          ! copy result back to b_cmplx
          b_cmplx = tmp_cmplx
    endif
endif


endif

end subroutine mul_coo_vec_cmplx
!==========================================================================================
  subroutine triplet_vector_cmplx_mult(A,x,b)
! multiply  A * x = b
  type (triplet), intent(in)               :: A
  complex(cwp),dimension(A%nj),intent(in)  :: x
  complex(cwp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_vector_cmplx_mult
  !==========================================================================================
  subroutine triplet_cmplx_vector_cmplx_mult(A,x,b)
! multiply  A * x = b
  type (triplet_cmplx), intent(in)               :: A
  complex(cwp),dimension(A%nj),intent(in)  :: x
  complex(cwp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_cmplx_vector_cmplx_mult
!==========================================================================================

!==========================================================================================
!---------------------------- TRIPLET INT MULTIPLY MATRIX----------------------------------
!==========================================================================================
subroutine mul_coo_int_mat(a_coo, m,n, b, method, c)
! C := A*b
implicit none

     type (triplet_int), intent(in)   :: a_coo
     integer, intent(in) :: m, n
  	 real(wp)				          :: b(m,n)
     real(wp)               :: c(a_coo%ni, n)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     real(wp):: alpha, beta


	if (a_coo%nj /= m) then
		print *, "Triplet(int)-matrix multiplications: dimensions must agree:", a_coo%nj, m
		stop
	end if

if (method == skit) then

!sparskit
         call triplet_int_matrix_mult(a_coo, m,n, b, c)
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

if (wp==dp) then
          call mkl_dcoomm(transa, a_coo%ni,n,m, alpha, matdescra, &
                          real(a_coo%vals, kind=dp), a_coo%indi, a_coo%indj, a_coo%nz, b, m, beta, c, a_coo%ni)
elseif (wp == sp) then
          call mkl_scoomm(transa, a_coo%ni,n,m, alpha, matdescra, &
                          real(a_coo%vals, kind=sp), a_coo%indi, a_coo%indj, a_coo%nz, b, m, beta, c, a_coo%ni)
endif


endif

end subroutine mul_coo_int_mat
!==========================================================================================
  subroutine triplet_int_matrix_mult(A, m, n, x,b)
! multiply  A * x = b
  type (triplet_int), intent(in)       :: A
  integer, intent(in) :: m, n
  real(wp),dimension(m,n), intent(in)      :: x
  real(wp),dimension(A%ni, n), intent(out) :: b
  integer :: k

	if (A%nj /= m) then
		print *, "Triplet(int)-matrix multiplications: dimensions must agree:", A%nj, m
		stop
	end if

  b = 0
  do k=1,n
    	call triplet_int_vector_mult(A,x(:,k),b(:,k))
  enddo

  end subroutine triplet_int_matrix_mult
!==========================================================================================
  subroutine triplet_int_vector_mult(A,x,b)
! multiply  A * x = b
  type (triplet_int), intent(in)       :: A
  real(wp),dimension(A%nj),intent(in)  :: x
  real(wp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_int_vector_mult
!==========================================================================================
!---------------------------- TRIPLET INT MULTIPLY MATRIX CMPLX----------------------------
!==========================================================================================
subroutine mul_coo_int_mat_cmplx(a_coo, m,n, b_cmplx, method, c_cmplx)
! C := A*b

implicit none

     type (triplet_int), intent(in)   :: a_coo
     integer, intent(in) :: m, n
  	 complex(cwp)				          :: b_cmplx(m,n)
     complex(cwp), allocatable            :: tmp_cmplx(:,:)
     complex(cwp)               :: c_cmplx(a_coo%ni, n)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     complex(cwp):: alpha, beta


	if (a_coo%nj /= m) then
		print *, "Triplet(int)-matrix multiplications: dimensions must agree:", a_coo%nj, m
		stop
	end if

if (method == skit) then

!sparskit
         call triplet_int_matrix_cmplx_mult(a_coo, m,n, b_cmplx, c_cmplx)
else

!mkl
transa = 'n'
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

if (cwp==dp) then
          call mkl_zcoomm(transa, a_coo%ni,n,m, alpha, matdescra, &
                          cmplx(a_coo%vals, 0, kind=dp), a_coo%indi, a_coo%indj, a_coo%nz, b_cmplx, m, beta, c_cmplx, a_coo%ni)
elseif (cwp == sp) then
          call mkl_ccoomm(transa, a_coo%ni,n,m, alpha, matdescra, &
                          cmplx(a_coo%vals, 0, kind=sp), a_coo%indi, a_coo%indj, a_coo%nz, b_cmplx, m, beta, c_cmplx, a_coo%ni)
endif

endif

end subroutine mul_coo_int_mat_cmplx
!==========================================================================================
  subroutine triplet_int_matrix_cmplx_mult(A, m, n, x,b)
! multiply  A * x = b
  type (triplet_int), intent(in)       :: A
  integer, intent(in) :: m, n
  complex(cwp),dimension(m,n), intent(in)      :: x
  complex(cwp),dimension(A%ni, n), intent(out) :: b
  integer :: k

	if (A%nj /= m) then
		print *, "Triplet(int)-matrix multiplications: dimensions must agree:", A%nj, m
		stop
	end if
  b = 0
  do k=1,n
    	call triplet_int_vector_cmplx_mult(A,x(:,k),b(:,k))
  enddo

  end subroutine triplet_int_matrix_cmplx_mult
!==========================================================================================
subroutine triplet_int_vector_cmplx_mult(A,x,b)
! multiply  A * x = b
  type (triplet_int), intent(in)       :: A
  complex(cwp),dimension(A%nj),intent(in)  :: x
  complex(cwp),dimension(A%ni),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
  	!print *, k, A%ni, A%nj, A%indi(k)
    b(A%indi(k)) = b(A%indi(k)) + A%vals(k) * x(A%indj(k))
  enddo

  end subroutine triplet_int_vector_cmplx_mult
!==========================================================================================
!---------------------------- MATRIX MULTIPLY TRIPLET INT ----------------------------------
!==========================================================================================
subroutine mul_mat_coo_int(m,n, b, a_coo, method, c)
! C := b*A

implicit none

     type (triplet_int), intent(in)   :: a_coo
     integer, intent(in) :: m, n
  	 real(wp), intent(in)	          :: b(m,n)
     real(wp), allocatable            :: b_tran(:,:), tmp(:,:)
     real(wp)               :: c(m, a_coo%nj)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     real(wp):: alpha, beta


	if (a_coo%ni /= n) then
		print *, "Triplet(int)-matrix multiplications: dimensions must agree:", a_coo%ni, n
		stop
	end if

if (method == skit) then

!sparskit
         call matrix_triplet_int_mult(m,n, b, a_coo, c)
else

!mkl
transa = 'T' ! transpose
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

allocate(b_tran(n,m), stat = istatus)
b_tran = transpose(b)
allocate(tmp(a_coo%nj,m), stat = istatus)

if (wp==dp) then
          call mkl_dcoomm(transa, a_coo%ni,m,a_coo%nj, alpha, matdescra, &
                          real(a_coo%vals, kind=dp), a_coo%indi, a_coo%indj, a_coo%nz, b_tran, n, beta, tmp, a_coo%nj)
elseif (wp == sp) then
          call mkl_scoomm(transa, a_coo%ni,m,a_coo%nj, alpha, matdescra, &
                          real(a_coo%vals, kind=sp), a_coo%indi, a_coo%indj, a_coo%nz, b_tran, n, beta, tmp, a_coo%nj)
endif
c = transpose(tmp)

endif

end subroutine mul_mat_coo_int
!==========================================================================================
  subroutine matrix_triplet_int_mult(m, n, x, A, b)
! multiply  A * x = b
  type (triplet_int), intent(in)       :: A
  integer, intent(in) :: m, n
  real(wp),dimension(m,n), intent(in)      :: x
  real(wp),dimension(m, A%nj), intent(out) :: b
  integer :: k

	if (A%ni /= n) then
		print *, "Matrix-triplet(int) multiplications: dimensions must agree:", A%ni, n
		stop
	end if

  b = 0
  do k=1,m
    call vector_triplet_int_mult(x(k,:), A, b(k,:))
  enddo

  end subroutine matrix_triplet_int_mult
!==========================================================================================
    subroutine vector_triplet_int_mult(x,A,b)
! multiply   x * A = b (x and b are row vectors)
  type (triplet_int), intent(in)       :: A
  real(wp),dimension(A%ni),intent(in)  :: x
  real(wp),dimension(A%nj),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indj(k)) = b(A%indj(k)) + A%vals(k) * x(A%indi(k))
  enddo

  end subroutine vector_triplet_int_mult
!==========================================================================================
!---------------------------- MATRIX CMPLX MULTIPLY TRIPLET INT ---------------------------
!==========================================================================================
subroutine mul_mat_cmplx_coo_int(m,n, b_cmplx, a_coo, method, c_cmplx)
! c_cmplx := b_cmplx*A
implicit none

     type (triplet_int), intent(in)   :: a_coo
     integer, intent(in) :: m, n
  	 complex(cwp), intent(in)	          :: b_cmplx(m,n)
     complex(cwp), allocatable            :: b_tran(:,:), tmp(:,:)
     complex(cwp)               :: c_cmplx(m, a_coo%nj)
     character :: method

     integer   :: istatus
     character        :: transa
     character(len=6) :: matdescra! = 'GOOF'
     complex(cwp):: alpha, beta

	if (a_coo%ni /= n) then
		print *, "Triplet(int)-matrix multiplications: dimensions must agree:", a_coo%ni, n
		stop
	end if

if (method == skit) then

!sparskit
		call matrix_cmplx_triplet_int_mult(m,n, b_cmplx, a_coo, c_cmplx)
else

!mkl
transa = 'T' ! transpose
alpha = 1.
beta = 0.
 matdescra(1:1) = 'G'
 matdescra(4:4) = 'F'

allocate(b_tran(n,m), stat = istatus)
b_tran = transpose(b_cmplx)
allocate(tmp(a_coo%nj,m), stat = istatus)

if (cwp==dp) then
          call mkl_zcoomm(transa,a_coo%ni,m,a_coo%nj, alpha, matdescra, &
                           cmplx(a_coo%vals, 0, kind=dp), a_coo%indi, a_coo%indj, a_coo%nz, b_tran, n, beta, tmp, a_coo%nj)
elseif (cwp == sp) then
          call mkl_ccoomm(transa, a_coo%ni,m,a_coo%nj, alpha, matdescra, &
                           cmplx(a_coo%vals, 0, kind=sp), a_coo%indi, a_coo%indj, a_coo%nz, b_tran, n, beta, tmp, a_coo%nj)
endif
c_cmplx = transpose(tmp)

endif

end subroutine mul_mat_cmplx_coo_int
!==========================================================================================

  subroutine matrix_cmplx_triplet_int_mult(m, n, x, A, b)
! multiply  A * x = b
  type (triplet_int), intent(in)       :: A
  integer, intent(in) :: m, n
  complex(cwp),dimension(m,n), intent(in)      :: x
  complex(cwp),dimension(m, A%nj), intent(out) :: b
  integer :: k

	if (A%ni /= n) then
		print *, "Matrix-triplet(int) multiplications: dimensions must agree:", A%ni, n
		stop
	end if

  b = 0
  do k=1,m
    call vector_cmplx_triplet_int_mult(x(k,:), A, b(k,:))
  enddo

  end subroutine matrix_cmplx_triplet_int_mult
!==========================================================================================
    subroutine vector_cmplx_triplet_int_mult(x,A,b)
! multiply   x * A = b (x and b are row vectors)
  type (triplet_int), intent(in)       :: A
  complex(cwp),dimension(A%ni),intent(in)  :: x
  complex(cwp),dimension(A%nj),intent(out) :: b
  integer :: k

  b = 0
  do k=1,A%nz
    b(A%indj(k)) = b(A%indj(k)) + A%vals(k) * x(A%indi(k))
  enddo

  end subroutine vector_cmplx_triplet_int_mult
!==========================================================================================

!---------------------------------- SORT TO COLUMN COMPLEX -----------------------------
!==========================================================================================
subroutine triplet_sort_col(Ait,Ajt,Ax,ipoint)
  real(wp),dimension(:),intent(inout)  :: Ax
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
!print *, "a"
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
!print *, "b"
  end subroutine triplet_sort_col
!==========================================================================================
!------------------------------------ SWAP COMPLEX------------------------------
!==========================================================================================
subroutine triplet_swap(k,l,Ait,Ajt,Ax,ipoint)
  real(wp),dimension(:),intent(inout)  :: Ax
  integer,dimension(:),intent(inout)            :: Ait,Ajt,ipoint
  integer,intent(in)          :: k,l
  integer                     :: itemp
  real(wp)           :: dtemp

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
!==========================================================================================
!-------------------- CONVERT TO Compressed Sparse Column ----------------------
!==========================================================================================
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
!==========================================================================================
!-------------------- CONVERT TO Compressed Sparse Column C-format -------------
!==========================================================================================
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
!==========================================================================================
!------------------------------------ SPYTRIPLET -------------------------------
!==========================================================================================
  subroutine triplet_spy(iwrite,Ait,Ajt,n)
  integer,dimension(:),intent(in)             :: Ait,Ajt
  integer,intent(in)                          :: iwrite
  character(len=1),dimension(:,:),allocatable :: elemnt
  character(len=1),dimension(:),allocatable   :: header
  integer :: istat,i,j,n,nz
  character(len=3) :: string

  nz   = size(Ait)

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
!==========================================================================================
!--------------------------- CSC1-MATRIX-SPY -----------------------------------
!==========================================================================================
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
!==========================================================================================
!==========================================================================================

end module my_sparse
