module my_lapack

     use precisions, only: sp, dp!, wp, cwp
	 use dispmodule
     use my_blas

!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!---------------------------- DOT PRODUCT ---------------------------------------
     interface gen_eig_nonsym
          module procedure GENERALIZED_EIGEN_NONSYM_sp
          module procedure GENERALIZED_EIGEN_NONSYM_dp
     end interface

     interface lapack_solve_general
          module procedure lapack_solve_general_sp
          module procedure lapack_solve_general_dp
     end interface  
!==========================================================================================

contains
!==========================================================================================
subroutine GENERALIZED_EIGEN_NONSYM_sp(A,B,n, cn, vn, eig_fun_flag)
! c := mat*vec
!or
! vec := mat*vec
	implicit none
!	include 'mkl_blas.fi'


     integer, intent(in)	:: N
	 real(sp)				:: A(N,N),B(N,N), cn(N)
     real(sp), allocatable	:: vn(:,:), vl(:,:)
     logical, optional		:: eig_fun_flag

     integer 				:: LDA,LDB,ldvl,ldvr, LWMAX, istat
     parameter				( LWMAX = 1000000 )
     real(sp), parameter	:: eps1 = 1e-12, eps2 = 1e-6
     integer :: ITYPE, INFO, LWORK
	 double precision, allocatable ::  WORK( : )
     real(sp)				:: alphar(N), alphai(N), beta(N)

     real :: START, FINISH
	 character*1	jobvl, jobvr

!     integer, optional :: n_threads
     integer :: num_threads = -1!, mkl_get_max_threads


!	num_threads = auto_num_threads
!	if (present(n_threads)) then
!		num_threads = n_threads
!	endif

	jobvl = 'N' ! the left generalized eigenvectors are not computed;
	ldvl = n
	jobvr = 'N' ! the right generalized eigenvectors are not computed by default
	ldvr = n
if (present(eig_fun_flag)) then
	if (eig_fun_flag) then
		jobvr = 'V' ! the right generalized eigenvectors are computed.
		ldvr = n
		allocate(vn(ldvr, n), stat=istat)
	endif
endif

!	matrix dims
        LDA=N
        LDB=N
!
    call CPU_TIME(START)

if (num_threads < 0) then
	! sequential lapack version

!     Query the optimal workspace.
        LWORK = -1
	allocate(work(1), stat=istat)
        call dggev( jobvl, jobvr, N, A, LDA, B, LDB, alphar, alphai, &
	                   beta, vl, ldvl, vn, ldvr, work, lwork, info)
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
	deallocate(work)

!     Solve eigenproblem.
	allocate(work(lwork), stat=istat)
        call DGGEV( jobvl, jobvr, N, A, LDA, B, LDB, alphar, alphai, &
	                   beta, vl, ldvl, vn, ldvr, work, lwork, info)
else
	! parallel version
!	if ( (num_threads >0).and.( num_threads /= mkl_get_max_threads()) ) then
!		call mkl_set_num_threads(num_threads) ! # threads
!	endif
endif

    call CPU_TIME(FINISH)

!     Check for convergence.
        if( INFO.GT.0 ) then
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        endif

!	calc the eigenvals
	    where (abs(alphai) > eps2) alphar = 0 ! if nonzero imag part in the eigenvals, discard
!	    	write(*,*) alphai
!	    	write(*,*) alphar

		if (count(abs(beta) < eps1) > 0) then
	    	write(*,*) 'Abs val < eps', (count(abs(beta) < eps1) > 0)
            write(*,*)'Computed generalized eigenvalues are unphysical. pnt 1'
            stop
        endif

	cn = alphar/beta

	    if ( count(cn > eps2) == 0 ) then
	    	write(*,*) 'Non-zero imag part', (vec_vec_dot(alphai,alphai) **.5 > eps2)
            write(*,*)'Computed generalized eigenvalues are unphysical. pnt 2'
            stop
        endif

end subroutine GENERALIZED_EIGEN_NONSYM_sp
!==========================================================================================
subroutine GENERALIZED_EIGEN_NONSYM_dp(A,B,n, cn, vn, eig_fun_flag)
! c := mat*vec
!or
! vec := mat*vec
	implicit none
!	include 'mkl_blas.fi'


     integer, intent(in)	:: N
	 real(dp)				:: A(N,N),B(N,N), cn(N)
     real(dp), allocatable	:: vn(:,:), vl(:,:)
     logical, optional		:: eig_fun_flag

     integer 				:: LDA,LDB,ldvl,ldvr, LWMAX, istat
     parameter				( LWMAX = 1000000 )
     real(dp), parameter	:: eps1 = 1e-12, eps2 = 1e-6
     integer :: ITYPE, INFO, LWORK
	 double precision, allocatable ::  WORK( : )
     real(dp)				:: alphar(N), alphai(N), beta(N)

     real :: START, FINISH
	 character*1	jobvl, jobvr

!     integer, optional :: n_threads
     integer :: num_threads = -1!, mkl_get_max_threads


!	num_threads = auto_num_threads
!	if (present(n_threads)) then
!		num_threads = n_threads
!	endif

	jobvl = 'N' ! the left generalized eigenvectors are not computed;
	ldvl = n
	jobvr = 'N' ! the right generalized eigenvectors are not computed by default
	ldvr = n
if (present(eig_fun_flag)) then
	if (eig_fun_flag) then
		jobvr = 'V' ! the right generalized eigenvectors are computed.
		ldvr = n
		allocate(vn(ldvr, n), stat=istat)
	endif
endif

!	matrix dims
        LDA=N
        LDB=N
!
    call CPU_TIME(START)

if (num_threads < 0) then
	! sequential lapack version

!     Query the optimal workspace.
        LWORK = -1
	allocate(work(1), stat=istat)
        call dggev( jobvl, jobvr, N, A, LDA, B, LDB, alphar, alphai, &
	                   beta, vl, ldvl, vn, ldvr, work, lwork, info)
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
	deallocate(work)

!     Solve eigenproblem.
	allocate(work(lwork), stat=istat)
        call DGGEV( jobvl, jobvr, N, A, LDA, B, LDB, alphar, alphai, &
	                   beta, vl, ldvl, vn, ldvr, work, lwork, info)
else
	! parallel version
!	if ( (num_threads >0).and.( num_threads /= mkl_get_max_threads()) ) then
!		call mkl_set_num_threads(num_threads) ! # threads
!	endif
endif

    call CPU_TIME(FINISH)

!     Check for convergence.
        if( INFO.GT.0 ) then
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        endif

!	calc the eigenvals
	    where (abs(alphai) > eps2) alphar = 0 ! if nonzero imag part in the eigenvals, discard
!	    	write(*,*) alphai
!	    	write(*,*) alphar

		if (count(abs(beta) < eps1) > 0) then
	    	write(*,*) 'Abs val < eps', (count(abs(beta) < eps1) > 0)
            write(*,*)'Computed generalized eigenvalues are unphysical. pnt 1'
            stop
        endif

	cn = alphar/beta

	    if ( count(cn > eps2) == 0 ) then
	    	write(*,*) 'Non-zero imag part', (vec_vec_dot(alphai,alphai) **.5 > eps2)
            write(*,*)'Computed generalized eigenvalues are unphysical. pnt 2'
            stop
        endif

end subroutine GENERALIZED_EIGEN_NONSYM_dp
!==========================================================================================

    ! =======
    ! Purpose
    ! =======
    !
    ! DGESV computes the solution to a real system of linear equations
    ! A * X = B,
    ! where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    !
    ! The LU decomposition with partial pivoting and row interchanges is
    ! used to factor A as
    ! A = P * L * U,
    ! where P is a permutation matrix, L is unit lower triangular, and U is
    ! upper triangular. The factored form of A is then used to solve the
    ! system of equations A * X = B.
    !
    ! Arguments
    ! =========
    !
    ! N (input) INTEGER
    ! The number of linear equations, i.e., the order of the
    ! matrix A. N >= 0.
    !
    ! NRHS (input) INTEGER
    ! The number of right hand sides, i.e., the number of columns
    ! of the matrix B. NRHS >= 0.
    !
    ! A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    ! On entry, the N-by-N coefficient matrix A.
    ! On exit, the factors L and U from the factorization
    ! A = P*L*U; the unit diagonal elements of L are not stored.
    !
    ! LDA (input) INTEGER
    ! The leading dimension of the array A. LDA >= max(1,N).
    !
    ! IPIV (output) INTEGER array, dimension (N)
    ! The pivot indices that define the permutation matrix P;
    ! row i of the matrix was interchanged with row IPIV(i).
    !
    ! B (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
    ! On entry, the N-by-NRHS matrix of right hand side matrix B.
    ! On exit, if INFO = 0, the N-by-NRHS solution matrix X.
    !
    ! LDB (input) INTEGER
    ! The leading dimension of the array B. LDB >= max(1,N).
    !
    ! INFO (output) INTEGER
    ! = 0: successful exit
    ! < 0: if INFO = -i, the i-th argument had an illegal value
    ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization
    ! has been completed, but the factor U is exactly
    ! singular, so the solution could not be computed.
!==========================================================================================
subroutine lapack_solve_general_sp(A,B)

    !solves Ax = b
    ! Can solve for multiple rhs. Change nrhs

    real(sp), intent(in) :: A(:,:)
    real(sp), intent(inout) :: B(:)

    integer :: n, i, j, info
    integer :: lda, ldb
    integer, allocatable :: ipiv(:)
    integer, parameter :: nrhs = 1

!    external :: dgesv

    n = size(A,1)
    lda = n
    ldb = n
    allocate(ipiv(n))

    CALL SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )


    ! Check for the exact singularity.
    IF( INFO.GT.0 ) THEN
	   WRITE(*,*)'The diagonal element of the triangular factor of A,'
       WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
       WRITE(*,*)'A is singular; the solution could not be computed.'
       stop
	elseif (info < 0) then
	   write(*,*)'Fejl i lapack_solve_general'
       write(*,*)'the i-th argument had an illegal value'
       stop
	endif

 end subroutine lapack_solve_general_sp
 !==========================================================================================
subroutine lapack_solve_general_dp(A,B)

    !solves Ax = b
    ! Can solve for multiple rhs. Change nrhs

    real(dp), intent(in) :: A(:,:)
    real(dp), intent(inout) :: B(:)

    integer :: n, i, j, info
    integer :: lda, ldb
    integer, allocatable :: ipiv(:)
    integer, parameter :: nrhs = 1

!    external :: dgesv

    n = size(A,1)
    lda = n
    ldb = n
    allocate(ipiv(n))

    CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )


    ! Check for the exact singularity.
    IF( INFO.GT.0 ) THEN
	   WRITE(*,*)'The diagonal element of the triangular factor of A,'
       WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
       WRITE(*,*)'A is singular; the solution could not be computed.'
       stop
	elseif (info < 0) then
	   write(*,*)'Fejl i lapack_solve_general'
       write(*,*)'the i-th argument had an illegal value'
       stop
	endif

 end subroutine lapack_solve_general_dp
!==========================================================================================
!==========================================================================================

end module my_lapack
