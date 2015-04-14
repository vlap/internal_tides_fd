module unsym_solvers

     use precisions, only: wp, cwp, dp, sp
     use iso_varying_string
     use my_sparse
     use dispmodule
     use control
!     use save_load, only: load_alloc_vec_int, save_vec_int


!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================

contains
!==========================================================================================

subroutine umfpack_unsym(n, csc, rhs, uvh, messages, filesave,loadsym, savenum)

!        use precisions, only: wp, cwp
!        use my_sparse

        implicit none

!       Input variables
     integer          :: n
     type (csc_cmplx) :: csc
     complex(cwp)     :: rhs(n)

     complex(cwp), allocatable     :: uvh(:)
     real(wp), allocatable         :: uvh_r(:), uvh_i(:)

     ! UMFPACK-variables
     real(wp), dimension(20) :: control
     real(wp), dimension(90) :: info
     integer                 :: numeric, symbolic, sys, status, filesave

     logical                 :: messages, savenum, loadsym

      ! set default parameters
        call umf4zdef (control)

      ! print control parameters.  set control (1) to 1 to print error messages only
      !                                set control (1) to 2 to print everything 
if (messages) then
      control(1) = 2
else
      control(1) = 1
endif

if (.not. loadsym) then

      call umf4zpcon (control)

      ! pre-order and symbolic analysis
      WRITE(*,'("Reordering ... ")', advance='no') 
      call umf4zsym(csc%ni, csc%nj, csc%indj, csc%indi, &
                   real(csc%vals, kind=wp), imag(csc%vals), symbolic, control, info)
!        call umf4zsym (n, n, Ap, Ai, Ax, Az, symbolic, control, info)
      WRITE(*,'("done. ")', advance='no')

!       print statistics computed so far
    if (messages) then
       call umf4zpinf (control, info)
    endif

!       check umf4zsym error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsym: ', info (1)
            stop
        endif

!       save the symbolic analysis to the file s1.umf
!       note that this is not needed until another matrix is
!       factorized, below.
        call umf4zssym (symbolic, filesave, status)
        if (status .lt. 0) then
            print *, 'Error occurred in umf4zssym: ', status
            stop
        endif

else
!       load the symbolic factorization back in (filename: n0.umf)
        call umf4zlsym (symbolic, filesave, status)
        if (status .lt. 0) then
            print *, 'Error occurred in umf4zlsym: ', status
            stop
        endif
endif

! numeric factorization
      WRITE(*,'("Factorization ... ")', advance='no') 
      call umf4znum(csc%indj, csc%indi, real(csc%vals, kind=wp), imag(csc%vals), &
                    symbolic, numeric, control, info)
!        call umf4znum (Ap, Ai, Ax, Az, symbolic, numeric, control, info)
      WRITE(*,'("done. ")', advance='no')

!       print statistics for the numeric factorization
    if (messages) then
       call umf4zpinf (control, info)
    endif

!       check umf4znum error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4znum: ', info (1)
            stop
        endif

if (savenum) then

!       save the LU factors to the file n0.umf
        call umf4zsnum (numeric, filesave, status)
        if (status .lt. 0) then
            print *, 'Error occurred in umf4zsnum: ', status
            stop
        endif
endif

!       free the symbolic analysis
        call umf4zfsym (symbolic)
!       free the numeric factorization
!        call umf4zfnum (numeric)
!       ----------------------------------------------------------------
!       in case to load the LU factors back in, and solve the system
!       ----------------------------------------------------------------
!       load the numeric factorization back in (filename: n0.umf)
!        call umf4zlnum (numeric, filesave, status)
!        if (status .lt. 0) then
!            print *, 'Error occurred in umf4zlnum: ', status
!            stop
!        endif

!       solve Ax=b, without iterative refinement
      WRITE(*,'("Solve ... ")', advance='no') 
      sys = 0
            allocate(uvh_r(n), uvh_i(n), stat=status)
      call umf4zsol(sys, uvh_r, uvh_i,  real(rhs, kind=wp), imag(rhs), numeric, control, info)
!        call umf4zsol (sys, x, xz, b, bz, numeric, control, info)
      WRITE(*,'("done. ")', advance='no')

        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
! solution
      allocate(uvh(n), stat=status)
            uvh = cmplx (uvh_r, uvh_i, kind=cwp)

!       free the numeric factorization
        call umf4zfnum (numeric)

!       No LU factors (symbolic or numeric) are in memory at this point.

!       print final statistics
if (messages) then
        call umf4zpinf (control, info)
endif


end subroutine umfpack_unsym

!==========================================================================================
!==========================================================================================

subroutine pardiso_unsym(n, csr, b, x, P, itn, CGS_ready, save_factors_dir) ! dir_save_descr

        !include 'mkl_pardiso.f90'
!        use precisions, only: wp, cwp, sp, dp

        implicit none

     type(params), intent(in) :: P

        !       Pardiso parameters
        integer, parameter :: mtype     = 13  ! complex unsymmetric
!        integer, parameter :: solver    =  0  ! use sparse direct method
        integer, parameter :: nrhs     = 1  ! Number of right-hand sides that need to be solved for.
        integer, parameter :: maxfct    = 1  ! Maximal number of factors with identical nonzero sparsity structure
                                             ! that the user would like to keep at the same time in memory.
                                             ! PARDISO can process several matrices with identical matrix sparsity pattern and
                                             ! store the factors of these matrices at the same time. Matrices with a different
                                             ! sparsity structure can be kept in memory with different memory address pointers pt.
        integer, parameter :: mnum    = 1  ! Actual matrix for the solution phase. With this scalar you can define the matrix
                                           !that you would like to factorize. The value must be: 1 ≤ mnum ≤ maxfct.

!       Internal solver memory pointer for 64-bit architectures
!       INTEGER*8 pt(64) Internal solver memory pointer (4 byte for 32-bit architectures, 8 byte for 64-bit architectures)
        integer*8, dimension(64)   :: pt = 0

!       Input variables
        type (csr_cmplx)           :: csr
        integer                    :: n, itn
!        logical					   :: last_itn_flag
        complex(cwp)               :: b(n)
        complex(cwp), allocatable  :: x(:)!, y(:)
        character(len=*),optional  :: save_factors_dir
!        character(len=*)		   :: dir_save_descr

        type(varying_string)	:: keep
        logical                 :: messages!, savenum, loadsym

!       All other variables
		logical					   :: CGS_ready
        integer                    :: phase, msglvl, istat, error
        integer                    :: iparm(64)
        integer                    :: idum !i
        integer                    :: perm
        complex(cwp)               :: ddum !, waltime1, waltime2

!		REAL*8					   :: dparm(64)
!									DPARM(33) Determinant for real symmetric indefinite matrices.
!									DPARM(34) Relative residual after Krylov-Subspace convergence.
!									DPARM(35) Number of Krylov-Subspace iterations.

select case (P%pardiso_symbolic)
	case (0)
		keep = 'none' !incurres symbolic factorization
	case (1)
		keep ='symbolic' !incurres symb+num factorization
	case default
		keep = 'none' !incurres symbolic factorization
end select

if (P%messages>=3) then
	messages = .true.
	msglvl = 1 ! print statistical information
else
	messages = .false.
	msglvl = 0 ! DO NOT print statistical information
endif

!  Setup Pardiso control parameters
	iparm = 0
	CALL pardiso_set_iparm(P, messages, CGS_ready, iparm)
	error = 0
!******** Direct, with parallel LU factorization *********
if ((iparm(4) == 0) .or. (.not. CGS_ready)) then

      if (((iparm(60)==1).or.(iparm(60)==2)) .and. (itn == 1)) then
			OPEN(unit=10, file='pardiso_ooc.cfg', action='write', status="replace") ! file=dir_save_descr//'/pardiso_ooc.cfg'
																					! have to export <MKL_PARDISO_OOC_CFG_PATH >/< MKL_PARDISO_OOC_CFG_FILE_NAME>
			if (present(save_factors_dir)) then
!				call system('export MKL_PARDISO_OOC_PATH=''' //  save_factors_dir // 'ooc_file''')
				write(10, '(a)') 'MKL_PARDISO_OOC_PATH = ' //  save_factors_dir // 'pardiso_ooc'
			endif
!			call system('export MKL_PARDISO_OOC_MAX_CORE_SIZE=' // tostring(P%pardiso_max_ram))
			write(10, '(a)') 'MKL_PARDISO_OOC_MAX_CORE_SIZE = ' // tostring(P%pardiso_max_ram)
			!call system('echo $MKL_PARDISO_OOC_MAX_CORE_SIZE=' // tostring(P%pardiso_max_ram))
!			call system('export MKL_PARDISO_OOC_MAX_SWAP_SIZE=' //  tostring(P%pardiso_max_swap))
			write(10, '(a)') 'MKL_PARDISO_OOC_MAX_SWAP_SIZE = ' //  tostring(P%pardiso_max_swap)
!			call system('export MKL_PARDISO_OOC_KEEP_FILE=1') ! 1 (default value), then all files are deleted, if it is set to 0, then all files are stored.
			write(10, '(a)') 'MKL_PARDISO_OOC_KEEP_FILE = ' //  tostring(P%pardiso_keep_files)

			close(10)
		endif

!   Do all steps (inc symbolic factorization) and free the memory

!   Initiliaze the internal solver memory pointer. Only necessary for the FIRST call of the PARDISO solver.
!   Reordering and Symbolic Factorization. This step also allocates all memory that is necessary for the factorization

  if ( (itn == 1).or.(trim(keep)=='none') ) then

! if both OOC mode and ITM modelling is on, then when switching Baro<->ITM, a new reordering does not clear the results of the old one.
! And things get messed up. As well as all available memory gets used up . Hence, an extra clean-up is needed before (potential) switching
	if (P%pardiso_ooc>=1) then
!		print *, "clearing up just in case"
		  phase     = -1           ! release ALL internal memory
	      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
	                    perm, nrhs, iparm, msglvl, ddum, ddum, error)
    endif

      phase     = 11      ! only reordering and symbolic factorization
      write(*, '(a)', advance='no') "   a) Reordering... "

    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, csr%vals, csr%indi, csr%indj, &
                    perm, nrhs, iparm, msglvl, ddum, ddum, error)
	call pardiso_error(error)

	if (P%messages >= 1) then
      write(*, '(a)') "done (peak-"//tostring(IPARM(15)/(1024.*1024.),'F0.1')// &
      			"Gb; permanent-"//tostring(IPARM(16)/(1024.*1024.),'F0.1')// "Gb) "
	endif
!      print *, " peak memory symbolic factorization ", iparm(15)/1024, "Mb"
!      print *, " permanent memory symbolic factorization ", iparm(16)/1024, "Mb"
!      print *, " estimate of size of factors ", iparm(17)/1024, "Mb"
!      print *, " Total peak memory consumption is ", real(max(IPARM(15), IPARM(16)+IPARM(17)))/(1024*1024), "Gb"

     if (messages) then
		print *, ""
     	   WRITE(*,*) 'Size of the factors   = ', iparm(17)/(1024.*1024.), "Gb"
           WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
           WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)
     endif

	if (P%messages >= 1) then
	!	Note: iparm(17), which is the size of the factors ONLY is used in the OCC criteria
	 	if (P%pardiso_ooc==1) then
	 	     if ( iparm(17)/1024 > (P%pardiso_max_ram + P%pardiso_max_swap) ) then
	 	     	if ( max(iparm(15), iparm(16) + iparm(63))/1024 < (P%pardiso_max_ram + P%pardiso_max_swap) ) then
	 	     		write(*, '(a)', advance='no') "      Size of matrix factors "
		    	 	WRITE(*,'(a)') tostring( iparm(17)/(1024.*1024.),'F0.1') // "Gb > "//&
	     		 	tostring((P%pardiso_max_ram + P%pardiso_max_swap)/1024.,'F0.1') //&
	     		 	"Gb RAM+Swap. Running out-of-core with peak memory usage: "//&
     			                tostring(max(iparm(15), iparm(16) + iparm(63))/(1024.*1024.),'F0.1') // "Gb. "
     		 	else
	     		 	WRITE(*,'(a)') "OCC memory usage: "//&
     			    tostring(max(iparm(15), iparm(16) + iparm(63))/(1024.*1024.),'F0.1') // "Gb. > "//&
					tostring((P%pardiso_max_ram + P%pardiso_max_swap)/1024.,'F0.1') //&
					"Gb RAM+Swap available. Not enough memory, matrix is too large."
					stop
				endif
     		 else
	 	     	write(*, '(a)', advance='no') "      Size of matrix factors "
	         	WRITE(*,'(a)') tostring(iparm(17)/(1024.*1024.),'F0.1') // "Gb < "//&
	     	 	tostring((P%pardiso_max_ram + P%pardiso_max_swap)/1024.,'F0.1') //"Gb RAM+Swap. "//&
	     	 	"Running in-core with peak memory usage: "//tostring(max(iparm(15),iparm(16)+iparm(17))/(1024.*1024.),'F0.1') // "Gb."
	    	endif
     	elseif (P%pardiso_ooc==0) then
	 	     if ( iparm(17)/1024 < (P%pardiso_max_ram + P%pardiso_max_swap) ) then
	 	     	 	write(*, '(a)', advance='no') "      Size of matrix factors "
		    	 	WRITE(*,'(a)') tostring( iparm(17)/(1024.*1024.),'F0.1') // "Gb > "//&
	     		 	tostring((P%pardiso_max_ram + P%pardiso_max_swap)/1024.,'F0.1')
	     			WRITE(*,'(a)') "      Running in-core. Peak memory usage: "//&
     			                tostring(max(iparm(15),iparm(16)+iparm(17))/(1024.*1024.),'F0.1') // "Gb."
	 	     else
     			WRITE(*,'(a)') "      Peak memory usage: "//&
     			                tostring(max(iparm(15),iparm(16)+iparm(17))/(1024.*1024.),'F0.1') // "Gb. > "//&
				tostring((P%pardiso_max_ram + P%pardiso_max_swap)/1024.,'F0.1') //"Gb RAM+Swap. Not enough memory, use out-of-core."
				stop
			endif
     	elseif (P%pardiso_ooc==2) then
	 	     if ( max(iparm(15), iparm(16) + iparm(63))/1024 < (P%pardiso_max_ram + P%pardiso_max_swap) ) then
     			WRITE(*,'(a)') "      Running out-of-core with peak memory usage: "//&
     			                tostring(max(iparm(15), iparm(16) + iparm(63))/(1024.*1024.),'F0.1') // "Gb. "//&
				"In-core version would require "//tostring(max(iparm(15),iparm(16)+iparm(17))/(1024.*1024.),'F0.1') // "Gb RAM+Swap."
	 	     else
     			WRITE(*,'(a)') "      OCC memory usage: "//&
     			                tostring(max(iparm(15), iparm(16) + iparm(63))/(1024.*1024.),'F0.1') // "Gb. > "//&
								tostring((P%pardiso_max_ram + P%pardiso_max_swap)/1024.,'F0.1') //&
								"Gb RAM+Swap available. Not enough memory, matrix is too large."
				stop
			 endif
		 endif
	endif
  endif

!   Factorization.
      phase     = 22  ! only factorization
      WRITE(*,'("   b) Factorization... ")', advance='no')
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, csr%vals, csr%indi, csr%indj, &
                    perm, nrhs, iparm, msglvl, ddum, ddum, error)
	if (P%messages >= 1) then
      WRITE(*,'(a)', advance='no') "done ("//tostring(iparm(17)/(1024.*1024.),'F0.1')//"Gb). "
	endif
	call pardiso_error(error)

!   Back substitution and iterative refinement
      phase     = 33  ! solve
!    Solve
      WRITE(*,'("   c) Solve... ")', advance='no')
      allocate(x(n), stat=istat)
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, csr%vals, csr%indi, csr%indj, &
                    perm, nrhs, iparm, msglvl, b, x, error)
	  call pardiso_error(error)
!	 WRITE(*,'(a, i2,a)', advance='no') "N of refinements ", iparm(7), ". "

!      WRITE(*,*) 'The solution of the system is '
!      DO i = 1, n
!        WRITE(*,*) ' x(',i,') = ', x(i)
!      END DO

!   Termination and release of memory
!if (P%pardiso_iterative == 0) then
  if (trim(keep)=='symbolic') then
      phase     = 0           ! release internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
                    perm, nrhs, iparm, msglvl, ddum, ddum, error)

  else!if (trim(keep)=='none') then
      phase     = -1           ! release ALL internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
                    perm, nrhs, iparm, msglvl, ddum, ddum, error)
  endif

! if OOC mode is on, then new initializing a new reordering does not clear the old one.
! And all available memory gets used up when switching Baro<->ITM. Hence, an extra RAM clean-up is needed before switching
!	if ( (phase==0) .and. ((P%pardiso_ooc==2).or.(P%pardiso_ooc==1).and.&
!					(max(iparm(15), iparm(16) + iparm(63))/1024 < (P%pardiso_max_ram + P%pardiso_max_swap))) ) then
!		print *, "last_itn_flag", last_itn_flag
!		if (last_itn_flag) then
!		  phase     = -1           ! release ALL internal memory
!	      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
!	                    perm, nrhs, iparm, msglvl, ddum, ddum, error)
!	    endif
!    endif

!endif
!******** End direct solver *********
else
!******** Iterative CGS (Conjugate-Gradients Squared) with LU factorization completed on step 1 **********
      phase     = 23  !   Iterations plus recomputation of factors L, U if necessary
      WRITE(*,'("   c) CGS iterations... ")', advance='no')
      allocate(x(n), stat=istat)
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, csr%vals, csr%indi, csr%indj, &
                    perm, nrhs, iparm, msglvl, b, x, error)

!	if (P%messages >= 1) then
!      WRITE(*,'(a)', advance='no') "done ("//tostring(iparm(17)/(1024.*1024.),'F0.1')//"Gb). "
!	endif
	if ((error .NE. 0).or.(iparm(20)<0)) then
		WRITE(*,*) "CGS failed, iparm(20) =", iparm(20)
	else
		! Report the number of iterations
		write (*,'(i2,a)', advance='no') iparm(20), ' iterations, '
    endif
    call pardiso_error(error)

!******** End iterative solver *********
endif

end subroutine pardiso_unsym

!==========================================================================================

!subroutine pardiso_clean(n, P, save_factors_dir)
!
!        implicit none
!
!        type(params) :: P
!
!        !       Pardiso parameters
!        integer, parameter :: mtype     = 13  ! complex unsymmetric
!        integer, parameter :: nrhs     = 1  ! Number of right-hand sides that need to be solved for.
!        integer, parameter :: maxfct    = 1  ! Maximal number of factors with identical nonzero sparsity structure
!                                             ! that the user would like to keep at the same time in memory.
!                                             ! PARDISO can process several matrices with identical matrix sparsity pattern and
!                                             ! store the factors of these matrices at the same time. Matrices with a different
!                                             ! sparsity structure can be kept in memory with different memory address pointers pt.
!        integer, parameter :: mnum    = 1  ! Actual matrix for the solution phase. With this scalar you can define the matrix
!                                           !that you would like to factorize. The value must be: 1 ≤ mnum ≤ maxfct.
!
       !Internal solver memory pointer for 64-bit architectures
       !INTEGER*8 pt(64) Internal solver memory pointer (4 byte for 32-bit architectures, 8 byte for 64-bit architectures)
!        integer*8, dimension(64)   :: pt = 0
!
       !Input variables
!        integer                    :: n
!        character(len=*),optional  :: save_factors_dir
       !All other variables
!        type(varying_string)	   :: keep
!        integer                    :: phase, error, msglvl, idum, perm
!        integer                    :: iparm(64) = 0
!        complex(cwp)               :: ddum
!
!
!select case (P%pardiso_symbolic)
!	case (0)
!		keep = 'none' !incurres symbolic factorization
!	case (1)
!		keep ='symbolic' !incurres symb+num factorization
!	case default
!		keep = 'none' !incurres symbolic factorization
!end select
!
  !Setup Pardiso control parameters
!	call pardiso_set_iparm(P, messages, CGS_ready, iparm)
!
   !Termination and release of memory
!	if (trim(keep)=='symbolic') then
!		phase     = -1           ! release ALL internal memory
!		CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
!                      perm, nrhs, iparm, msglvl, ddum, ddum, error)
!	endif
!
!	if ((P%pardiso_ooc>=1).and.(present(save_factors_dir))) then
!		call system('find ' // save_factors_dir // ' -iname ' // '''pardiso_ooc.*''' // ' -exec rm -f {} \;' )
!	endif
!
!end subroutine pardiso_clean

!==========================================================================================
subroutine pardiso_set_iparm(P, messages, CGS_ready, iparm)

        implicit none

        type(params), intent(in) :: P
        logical, intent(in)      :: messages, CGS_ready
        integer      			 :: iparm(64)
		integer                  :: iterative, OMP_GET_MAX_THREADS

iparm(64) = 0

!	INPUT PARAMETERS
      iparm(1) = 1 ! no solver default values. If iparm(1) = 0, iparm(2) through iparm(64) are filled with default values
		if ((P%omp_num_threads > 0).and.(P%omp_num_threads < OMP_GET_MAX_THREADS())) then
			iparm(3) = P%omp_num_threads ! numbers of processors, value of OMP_NUM_THREADS
	    else
			iparm(3) = OMP_GET_MAX_THREADS() ! numbers of processors, value of OMP_NUM_THREADS
	    endif
	    if (iparm(3) > 1) then
			iparm(2) = 3   ! fill-in reducing ordering: 0, the minimum degree algorithm is applied
	        	           ! 2, the solver uses the nested dissection algorithm from the METIS package
	            	       ! 3, the parallel (OpenMP) version of the nested dissection algorithm is used.
	            	       !    It can decrease the time of computations on multi-core computers, especially when PARDISO Phase 1 takes significant time.
		else
			iparm(2) = 2
		endif
      iparm(5) = 0 ! 0 -- no user fill-in reducing permutation, save computed perm vector for future iterations
                   ! 1 -- user supplies the permutation vector perm
                   ! 2 -- save computed permutation vector into perm
      iparm(6) = 0 ! 0 solution on the first n compoments of x
      iparm(8) = 2 ! numbers of iterative refinement steps
      			   ! If iparm(8)< 0, the accumulation of the residue uses extended precision real and complex data types.
				   ! Perturbed pivots result in iterative refinement (independent of iparm(8)=0) and the number of executed iterations is reported in iparm(7).
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13

!************ FOLLOWING TWO PARAMS SHOULD BOTH BE 0 FOR FASTER COMPUTATIONS ***************************
!***************** with LESS RAM CONSUMPTION ON THE REORDERING STAGE ***************************
! Must provide the numerical values of the matrix in the analysis phase (PHASE=11)
! if IPARM(11) (scaling) or PARM(13)(matchings) are NOT zeros.

      iparm(11) = 0 ! Scaling for nonsymmetric matrices (MTYPE=11 or MTYPE=13): permute large elements close the diagonal and to scale the matrix so that the largest elements are equal to 1.
      				! 0 -- do not use nonsymmetric permutation and scaling MPS (faster)
                    ! 1 -- use nonsymmetric permutation and scaling MPS (more accurate and stable)
      iparm(13) = 0 ! Improved accuracy using (non-)symmetric weighted matchings.
      				! 0 -- maximum weighted matching algorithm is switched-off (faster)
                    ! 1 -- maximum weighted matching algorithm is switched-on (more accurate and stable)
                    ! 2 -- advanced matchings, higher accuracy
!*******************************************************************************************************
      iparm(12) = 0 ! solving with transposed or conjugate transposed matrix (1 -- solve with conjugate transposed A, 2 -- solve with transposed A)

	  if (iparm(3) >= 8) then
	      iparm(24) = 1 ! parallel factorization control: 1 -- new two-level scheduling algorithm. this algorithm generally improves scalability in case of parallel factorization on many threads (>= 8).
	  else
	      iparm(24) = 0 ! parallel factorization control: PARDISO uses the previous parallel factorization.
	  endif
	  if (iparm(3) > 1) then
	      iparm(25) = 0 ! parallel forward/backward solve control: PARDISO uses a parallel algorithm for the solve step.
	  else
	      iparm(25) = 1 ! PARDISO uses sequential forward and backward solve. (IC mode only)
	  endif
	  iparm(26) = 0	! Splitting of Forward/Backward Solve: 0 indicates that a normal solve step is performed
	  iparm(27) = 0 ! 1 -- matrix consistency checker (whether column indices are sorted in increasing order within each row.)
	! single or double precision of PARDISO: can be changed only during the solver's phase 1.
      if (cwp==dp) then
           iparm(28) = 0 !for double precision pardiso
      elseif (cwp==sp) then
           iparm(28) = 1 !for single precision pardiso
	  else
	  		print *, "Error in pardiso_unsym (precisions)"
	  		stop
      endif
      !iparm(29) = 0 ! Switch between 32-bit and 64-bit factorization.
      !iparm(30) = 0 ! Control the size of the supernodes.
      !iparm(31) = 0 ! partial solution for sparse right-hand sides and sparse solution: 0 (default value), this option is disabled.
      !iparm(32) = 0 ! Use the multi-recursive iterative linear solver.
      !iparm(33) = 0 ! Determinant of a real symmetric indefinite matrix.
      !iparm(34) = 0 ! Identical solution independent on the number of processors.
      !iparm(35) = 0 ! C or Fortran style array indexing.
      !IPARM(35) to IPARM(50) | Unused.
!	Parallel Distributed-Memory options
	  if (P%mpi_num_nodes > 1) then
	  	iparm(51) = 1
	  	iparm(52) = P%mpi_num_nodes
	  endif
      !IPARM(53) to IPARM(60) | Unused.

!************* SWITCH BETWEEN DIRECT AND ITERATIVE SOLVER **************
if ((P%pardiso_iterative == 1).and.(CGS_ready)) then
	iterative = 1 ! iterative solver is chosen, LU preconditioner for the CGS iterative solver is calculated on itn=1
else
	iterative = 0 ! direct solver
endif

	select case (iterative)
		case (0) ! direct solver
		  iparm(4) = 0 ! no iterative algorithm (preconditioned CGS)
			! Out-of-core (OC) version or in-core (IC) version. The OC PARDISO can solve very large problems by holding the matrix factors in files on the disk (the amount of main memory required is reduced).
	      if ((P%pardiso_ooc==1).or.(P%pardiso_ooc==2)) then
	           iparm(60) = P%pardiso_ooc ! IC PARDISO is used if iparm(60) = 1 and the total memory of RAM (in Mb) is insufficient. If iparm(60) = 2 then OCC mode is used.
	!-----------TURN OFF FEATURES NOT AVAILABLE IN OUT-OF-CORE MODE------------!
				iparm(5) = 0  ! (user permutation)
				iparm(8) = 0  ! (iterative refinement steps)
				iparm(24) = 0 ! (two-level factorization algorithm)
				iparm(31) = 0 ! (partial solution for sparse right-hand sides and sparse solution)
	      endif

		case (1) ! iterative solver (preconditioned CGS)
		! iparm(4) has the form iparm(4)= 10*L+K.
		! K=1:	CGS iteration replaces the computation of LU.
		! The preconditioner is LU that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
		! The value L controls the stopping criterion of the Krylow-Subspace iteration: epsCGS = 10^{-L} is used in the stopping criterion
           iparm(4) = 10*P%pardiso_iter_tol + 1 !LU-preconditioned CGS iteration with a stopping tolerance of 10^{-L} (1 for nonsymmetric matrices)
           iparm(60) = 0 ! Turn off OCC, otherwise error pops up
	    case default
	        print *, "Error in pardiso_unsym (SWITCH BETWEEN DIRECT AND ITERATIVE SOLVER)"
	  		stop
	 end select

!	OUTPUT PARAMETERS
      iparm(7) = 0  ! number of performed iterative refinement steps that are actually performed during the solve step.
      iparm(14) = 0 ! Output: number of perturbed pivots during the elimination process for mtype =11, mtype =13, mtype =-2, mtype =-4, or mtype =-6.
      iparm(15) = 0 ! peak memory symbolic factorization in (in Kb, computed in phase 1)
      iparm(16) = 0 ! permanent memory symbolic factorization (in Kb, computed in phase 1)
      iparm(17) = 0 ! size of factors / memory numerical factorization and solution (in Kb, computed in phase 1, IC mode only)
      if (messages) then
	      iparm(18) = -1 ! number of nonzeros in the factor LU (if iparm(18) < 0 on entry)
	      iparm(19) = -1 ! number of operations in MFLOPS that are necessary to factor the matrix A (if iparm(19) < 0 on entry)
      else
	      iparm(18) = 0 ! number of nonzeros in the factor LU (if iparm(18) < 0 on entry)
	      iparm(19) = 0 ! number of operations in MFLOPS that are necessary to factor the matrix A (if iparm(19) < 0 on entry)
      endif
      iparm(20) = 0 ! CG/CGS diagnostics. If iparm(20)> 0, iparm(20) = number of CG iterations, otherwise an error code
	  iparm(22) = 0 ! number of positive eigenvalues for SYMMETRIC indefinite matrices.
	  iparm(23) = 0 ! number of negative eigenvalues for SYMMETRIC indefinite matrices.
	  iparm(30) = 0 ! number of the equation where PARDISO detects zero or negative pivot (for matrix types mtype = 2 and mtype = 4)
	  iparm(63) = 0 ! size of the minimum OOC memory for numerical factorization and solution (in Kb, computed in phase 1, OCC mode only)

end subroutine pardiso_set_iparm
!==========================================================================================
subroutine pardiso_error(error)

        implicit none
		integer                  :: error

select case (error)
	case (0)
		return;
	case (-1)
		WRITE(*,*) 'ERROR: '// 'Input inconsistent.'
	case (-2)
		WRITE(*,*) 'ERROR: '// 'Not enough memory. '//&
		            'Consider changing iparm(11) and iparm(13) to 0 to reduce RAM consumption on reordering stage.'
	case (-3)
		WRITE(*,*) 'ERROR: '// 'Reordering problem.'
	case (-4)
		WRITE(*,*) 'ERROR: '// 'Zero pivot, numerical fact. or iterative refinement problem.'
	case (-5)
		WRITE(*,*) 'ERROR: '// 'Unclassified (internal) error.'
	case (-6)
		WRITE(*,*) 'ERROR: '// 'Preordering failed (matrix types 11, 13 only).'
	case (-7)
		WRITE(*,*) 'ERROR: '// 'Diagonal matrix problem.'
	case (-8)
		WRITE(*,*) 'ERROR: '// '32-bit integer overflow'
	case (-9)
		WRITE(*,*) 'ERROR: '// 'not enough memory for OOC'
	case (-10)
		WRITE(*,*) 'ERROR: '// 'No license file pardiso.lic found.'
	case (-11)
		WRITE(*,*) 'ERROR: '// 'License is expired.'
	case (-12)
		WRITE(*,*) 'ERROR: '// 'Wrong username or hostname.'
	case (-100)
		WRITE(*,*) 'ERROR: '// 'Reached maximum number of Krylov-subspace iteration in iterative solver.'
	case (-101)
		WRITE(*,*) 'ERROR: '// 'No sufficient convergence in Krylov-subspace iteration within 25 iterations.'
	case (-102)
		WRITE(*,*) 'ERROR: '// 'Error in Krylov-subspace iteration.'
	case (-103)
		WRITE(*,*) 'ERROR: '// 'Break-Down in Krylov-subspace iteration.'
!********************************************************************
	case default
		WRITE(*,*) 'Unknown ERROR was detected: ', error
end select

IF (error .NE. 0) THEN
!        STOP 1
END IF

end subroutine pardiso_error
!==========================================================================================

end module unsym_solvers
