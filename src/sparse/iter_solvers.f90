module iter_solvers

     use precisions, only: wp, cwp, dp, sp
!     use iso_varying_string
     use my_sparse
     use my_sparse_aggregate
     use dispmodule
     use control
     use save_load, only: !save_sparse


!==========================================================================================
!**************************** interfaces **************************************************
!==========================================================================================

contains
!==========================================================================================

subroutine mkl_gmres(n, csr, b, x, P, itn, cpt, dir_sols)

!************** WARNING: the routine dfgmres works in DOUBLE PRECISION only ***************
        !include 'mkl_pardiso.f90'
!        use precisions, only: wp, cwp, sp, dp

      implicit none

      include "mkl_rci.fi"

      integer n, size
      parameter (size=128)			! parameter arrays are of size 128.

!       Input variables
        type (csr_cmplx)           :: csr
        complex(cwp), allocatable  :: b(:)!, expected_solution(:)
        integer                    :: itn
        character(len=2) :: cpt
        character(len = *) :: dir_sols

		type(params), intent(in)   :: P
!---------------------------------------------------------------------------
! Define arrays for a (2n x 2n) system of eqs with REAL matrix and rhs
! Compressed sparse row storage is used for sparse representation
! Allocate storage for the ipar parameters and the solution/rhs/residual vectors
!---------------------------------------------------------------------------
        type(csr_dp)		       :: a_re, prec_ilut, &	! a_re - real-valued analog of the system csr*x = b, then
        								prec_ilutp		! a_re is a (2n x 2n) matrix, a_re = [Re(csr) -Im(csr); Im(csr) Re(csr)]
		type(triplet_dp)		   :: a_re_coo
        real(dp), allocatable  	   :: b_re_copy(:), residual(:)
        integer			:: filled, nzmax

        real(dp), allocatable  	   :: b_re(:), trvec(:)
        integer                    :: ipar(size),ierr
        real(dp)			  	   :: dpar(size)
        real(dp)			  	   :: tmp((2*n)*(2*(2*n)+1)+((2*n)*((2*n)+9))/2+1)
!        real(dp)			  	   :: tmp(((2*P%gmres_rest_it+1)*n+P%gmres_rest_it*P%gmres_rest_it+9)/2 + 1)) ! array of size ((2*ipar(15)+1)*n+ipar(15)*ipar(15)+9)/2 + 1), ipar(15)=P%gmres_rest_it

        real(dp), allocatable  	   :: jw(:), w(:)
        integer, allocatable  	   :: iperm(:)

        integer                    :: incx
		parameter (incx=1)
!        real(dp)			  	   :: nrm2


        real(dp), allocatable		:: computed_solution(:)
        complex(cwp), allocatable	:: u(:), v(:), h(:)
        integer				 		:: nu, nv, nh
!       Output
        complex(cwp), allocatable  :: x(:)
!---------------------------------------------------------------------------
! some additional variables to use with the rci (p)fgmres solver
!---------------------------------------------------------------------------
		integer                    :: itercount, rci_request, j, istat
        real(dp)			  	   :: dvar
!        logical					   :: messages
!---------------------------------------------------------------------------
! some additional variables to use for ilut preconditioner call
!---------------------------------------------------------------------------
		integer                    :: maxfil
        real(dp)			  	   :: tol, permtol
!---------------------------------------------------------------------------
! an external blas function is taken from mkl blas to use
! with the rci (p)fgmres solver
!---------------------------------------------------------------------------
      double precision dnrm2, dznrm2
      external dnrm2, dznrm2

! Compare the result with PARDISO
!		allocate(expected_solution(n), stat = istat)
!		expected_solution = x
!		deallocate(x)

	write(*, '(a)', advance='no') "   a) Initializing "
!---------------------------------------------------------------------------
! initialize the real-valued matrix a_re and right hand side b_re
!---------------------------------------------------------------------------
		! Matrix a_re is a (2n x 2n) matrix, a_re = [Im(csr) Re(csr); Re(csr) -Im(csr)]
		call conv_csr_cmplx2csr_dp(csr, a_re)
!		call csr2coo(a_re, a_re_coo)
!	      call save_sparse(a_re_coo, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/'//&
!	      							  'baro_fd/out/0000_00_00__00_00/global/mats/gmres.dat')
!	      print *, "Saved a_re_coo with nnz", a_re_coo%nz

!		now the rhs
		allocate(b_re(2*n), stat = istat)
		b_re = (/ real(imag(b),dp), real(b, dp)/) ! b_re = [Im(b) Re(b)]
		deallocate(b)
!---------------------------------------------------------------------------
! save the right-hand side b_re in vector b_re_copy for future use
!---------------------------------------------------------------------------
		allocate(b_re_copy(2*n), stat = istat)
		call dcopy(2*n, b_re, 1, b_re_copy, 1)
!---------------------------------------------------------------------------
! initialize the initial guess. use the solution computed on the previous step
!---------------------------------------------------------------------------
		allocate(computed_solution(2*n), stat = istat)

		if (itn>1) then
	     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat', nu)
	     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat', nv)
	     call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat', nh)
	     computed_solution = (/real(u, dp),real(v, dp),real(h, dp), real(imag(u),dp), real(imag(v),dp),real(imag(h),dp)/)
	     deallocate(u,v,h)
		else
	     computed_solution=0.d0
	    endif

!---------------------------------------------------------------------------
! initialize the solver
!---------------------------------------------------------------------------
      call dfgmres_init(2*n, computed_solution, b_re, rci_request, ipar, dpar, tmp)
      if (rci_request.ne.0) goto 999
!---------------------------------------------------------------------------

	write(*, '(a)', advance='no') "and calculating "//trim(P%gmres_prec)//" preconditioner..."
!---------------------------------------------------------------------------
! calculate ilut preconditioner.
!                      !attention!
! dcsrilut routine uses some ipar, dpar set by dfgmres_init routine.
! important for dcsrilut default entries set by dfgmres_init are
!
! if ilut is going to be used out of mkl fgmres context, than the values
! of ipar(2), ipar(6), ipar(31), and dpar(31), should be user
! provided before the dcsrilut routine call.
!---------------------------------------------------------------------------
	itercount = 1

	if (P%messages>=3) then
		ipar(2) = 6	! output of ALL error and warning messages to the screen,
	elseif (P%messages>=1) then
		ipar(2) = 1	! the error and warning messages are written to the newly created file MKL_RCI_FGMRES_Log.txt
	else
		ipar(6) = 0	! error and warning messages are not generated at all.
		ipar(7) = 0 ! error and warning messages are not generated at all.
	endif

! ipar(31)= 0 - abort dcsrilut calculations if routine meets zero diagonal element.
	ipar(31)=1		! ipar(31)= 1 - change small diagonal value to that given by dpar(31)
	dpar(31)=1.d-5	! dpar(31)= 1.d-5  instead of the default value set by dfgmres_init.
!                   it is the target value of the diagonal value if it is
!                   small as compared to given tolerance multiplied
!                   by the matrix row norm and the routine should
!                   change it rather than abort dcsrilut calculations.
	tol=1.d-6		! Tolerance for threshold criterion for the resulting entries of the preconditioner.
	maxfil=n	! Maximum fill-in, which is half of the preconditioner bandwidth.
					! The number of non-zero elements in the rows of the preconditioner can not exceed (2*maxfil+1).
	nzmax=(2*maxfil+1)*a_re%ni - maxfil*(maxfil+1)+1 ! Max number of non-zeros in the preconditioner
							! - maxfil*(maxfil+1) -- number of elements cut off from the banded matrix by the "corners"
!---------------------------------------------------------------------------
    call init_sparse_0(prec_ilut, a_re%ni,a_re%nj, nzmax)

	call dcsrilut(2*n, a_re%vals, a_re%indi, a_re%indj, prec_ilut%vals, prec_ilut%indi, prec_ilut%indj, tol, maxfil, ipar, dpar, ierr)
!	nrm2=dnrm2(prec_ilut%nz, prec_ilut%vals, incx)

	if(ierr.ne.0) then
	  write(*,'(a,i4)') ' error after calculation of the preconditioner dcsrilut',ierr
	  goto 998
	elseif (P%messages >= 1) then
      write(*, '(a)') " done"
	endif
!---------------------------------------------------------------------------
!-------------------------	SPARSEKIT ILUTP --------------------------------
!---------------------------------------------------------------------------
if (1==0) then
	tol=1.d-6		! Tolerance for threshold criterion for the resulting entries of the preconditioner.
	maxfil = n !20
	permtol = 1	! tolerance ratio used to  determne whether or not to permute two columns.
					!  At step i columns i and j are permuted when abs(a(i,j))*permtol .gt. abs(a(i,i))
					! [0 --> never permute; good values 0.1 to 0.01]
!=============
! work arrays:
!=============
    call init_sparse_0(prec_ilutp, a_re%ni,a_re%nj, (2*maxfil+1)*a_re%ni - maxfil*(maxfil+1)+1 )

	allocate(jw(4*n), stat = istat)		! jw      = integer work array of length 2*N (N=2n here)
	allocate(w(2*n+1), stat = istat)	! w       = real work array of length N+1 (N=2n here)
	allocate(iperm(4*n), stat = istat)	!contains the permutation arrays.
										!           iperm(1:n) = old numbers of unknowns
										!           iperm(n+1:2*n) = reverse permutation = new unknowns.

	call ilutp(2*n, a_re%vals, a_re%indi, a_re%indj,maxfil,tol,permtol,&
			   2*n, prec_ilutp%vals, prec_ilutp%indi, prec_ilutp%indj, nzmax,w,jw,iperm,ierr)
	print *, ""
	print *, "hi"
	print *, ""
endif

!---------------------------------------------------------------------------
! set the desired parameters:
! do the restart after P%gmres_rest_it iterations
      ipar(15)=P%gmres_rest_it
! do not do the stopping test for the maximal number of iterations
      ipar(8)=0
! do the preconditioned iterations of fgmres method
      ipar(11)=1
! set the relative tolerance to 10^{-P%gmres_tol}
      dpar(1)=10**(-P%gmres_tol)
!---------------------------------------------------------------------------
! check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
      call dfgmres_check(2*n, computed_solution, b_re, rci_request, ipar, dpar, tmp)
      if (rci_request.ne.0) goto 999
!---------------------------------------------------------------------------
! print the info about the rci fgmres method
!---------------------------------------------------------------------------
	if (P%messages>=3) then
		call mkl_gmres_print_info(ipar)
	endif

	if (P%messages<3) then
		WRITE(*,'("   c) Solve... ")', advance='no')
    else
    	WRITE(*,'("   c) Solve... ")')
   	endif
!---------------------------------------------------------------------------
! compute the solution by rci (p)fgmres solver with preconditioning
! reverse communication starts here
!---------------------------------------------------------------------------
 1     call dfgmres(2*n, computed_solution, b_re_copy, rci_request, ipar, dpar, tmp)
! 		call disp("1")
!---------------------------------------------------------------------------
! if rci_request=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
      if (rci_request.eq.0) goto 3
!---------------------------------------------------------------------------
! if rci_request=1, then compute the vector a*tmp(ipar(22))
! and put the result in vector tmp(ipar(23))
!---------------------------------------------------------------------------
      if (rci_request.eq.1) then
      	call mkl_dcsrgemv('n',a_re%ni, a_re%vals, a_re%indi, a_re%indj, tmp(ipar(22)), tmp(ipar(23)))
      	goto 1
      endif
!---------------------------------------------------------------------------
! if rci_request=2, then do the user-defined stopping test
! the residual stopping test for the computed solution is performed here
!---------------------------------------------------------------------------
! note: from this point vector b_re(n) is no longer containing the right-hand
! side of the problem! it contains the current fgmres approximation to the
! solution. The right-hand side is saved in b_re_copy
!---------------------------------------------------------------------------
      if (rci_request.eq.2) then
! request to the dfgmres_get routine to put the solution into b_re
      	ipar(13)=1	! the routine writes the solution to the right hand side vector b_re
! get the current fgmres solution in the vector b_re
      	call dfgmres_get(2*n, computed_solution, b_re, rci_request, ipar, dpar, tmp, itercount)
! compute the current true residual via mkl (sparse) blas routines
		if (.not. allocated(residual)) allocate(residual(2*n), stat = istat)
      	call mkl_dcsrgemv('n', a_re%ni, a_re%vals, a_re%indi, a_re%indj, b_re, residual)
      	call daxpy(2*n, -1.0d0, b_re_copy, 1, residual, 1)
      	dvar=dnrm2(2*n, residual, 1)

		if (P%messages>=3) then
			WRITE(*,'(a)', advance='no') "itn "//tostring(itercount)//", ||res||="//tostring(dvar)//"; "
		endif

      	if (dvar .lt. 10.**(-P%gmres_tol)) then
!	   		print *, dvar, " < ", 10.**(-P%gmres_tol)
      	   goto 3
      	else
      	   goto 1
      	endif
      endif
!---------------------------------------------------------------------------
! if rci_request=3, then apply the preconditioner on the vector
! tmp(ipar(22)) and put the result in vector tmp(ipar(23))
! here is the recommended usage of the result produced by ilut routine
! via standard mkl sparse blas solver routine mkl_dcsrtrsv.
!---------------------------------------------------------------------------
      if (rci_request.eq.3) then
!		if (P%messages>=3) then
!			WRITE(*,'(a)', advance='no') "itn "//tostring(itercount)//", applying precond; "
!		endif
      	if (.not. allocated(trvec)) allocate(trvec(2*n), stat = istat)
       	call mkl_dcsrtrsv('l','n','u',2*n,prec_ilut%vals,prec_ilut%indi,prec_ilut%indj, tmp(ipar(22)),trvec)
       	call mkl_dcsrtrsv('u','n','n',2*n,prec_ilut%vals,prec_ilut%indi,prec_ilut%indj, trvec,tmp(ipar(23)))
       	goto 1
      endif
!---------------------------------------------------------------------------
! if rci_request=4, then check if the norm of the next generated vector is
! not zero up to rounding and computational errors. the norm is contained
! in dpar(7) parameter
!---------------------------------------------------------------------------
!	Assuming that 1.0d-16 is as far as we can/want to go
      if (rci_request.eq.4) then
      	if (dpar(7).lt.1.0d-16) then
      	   goto 3
      	else
      	   goto 1
      	endif
!---------------------------------------------------------------------------
! if rci_request=anything else, then dfgmres subroutine failed
! to compute the solution vector: computed_solution(n)
!---------------------------------------------------------------------------
      else
      	goto 999
      endif
!---------------------------------------------------------------------------
! reverse communication ends here
! get the current iteration number and the fgmres solution. (do not forget to
! call dfgmres_get routine as computed_solution is still containing
! the initial guess!). request to dfgmres_get to put the solution into
! vector computed_solution(n) via ipar(13)
!---------------------------------------------------------------------------
 3     ipar(13)=0
      call dfgmres_get(n, computed_solution, b_re_copy, rci_request, ipar, dpar, tmp, itercount)

      allocate(x(n), stat = istat)
      x = cmplx(computed_solution(1:n), computed_solution(n+1:2*n), kind=cwp)
!---------------------------------------------------------------------------
! print solution vector: computed_solution(n) and
! the number of iterations: itercount
!---------------------------------------------------------------------------
	if (P%messages >= 1) then
		write( *,'(a,i2,a)', advance='no') 'done after ',itercount, ' iterations '
	endif
!      write( *,'(a)', advance='no') 'L2 norm of the difference between the expected and the obtained solutions:'
!      write(*,'(3e14.7)') dznrm2(n, x - expected_solution, 1)/dznrm2(n, expected_solution, 1)
      goto 1000

999	if (P%messages >= 1) then
		write( *,'(a,i2)') 'the solver has returned the error code ', rci_request
	endif
998 if (P%messages >= 1) then
		write( *,'(a,a)') 'unfortunately, fgmres+ilut fortran example has failed'
	endif

1000  continue

end subroutine mkl_gmres

!==========================================================================================

subroutine mkl_gmres_print_info(ipar)

      implicit none

        integer, intent(in) :: ipar(:)

      write( *,'(a)') ' '
      write( *,'(a,a)') 'some info about the current run of rci fgmres', &
       ' method:'
      write( *,'(a)') ' '
      if (ipar(8).ne.0) then
         write(*,'(a,i1,a,a)') 'as ipar(8)=',ipar(8),', the automatic', &
       ' test for the maximal number of iterations will be performed'
      else
      	write(*,'(a,i1,a,a)') 'as ipar(8)=',ipar(8),', the automatic', &
         ' test for the maximal number of iterations will be skipped'
      endif
      write( *,'(a)') '+++'
      if (ipar(9).ne.0) then
      	write(*,'(a,i1,a,a)') 'as ipar(9)=',ipar(9),', the automatic', &
       ' residual test will be performed'
      else
      	write(*,'(a,i1,a,a)') 'as ipar(9)=',ipar(9),', the automatic', &
       ' residual test will be skipped'
      endif
      write( *,'(a)') '+++'
      if (ipar(10).ne.0) then
      	write(*,'(a,i1,a,a)') 'as ipar(10)=',ipar(10),', the', &
       ' user-defined stopping test will be requested via rci_request=2'
      else
      	write(*,'(a,i1,a,a,a)') 'as ipar(10)=',ipar(10),', the', &
       ' user-defined stopping test will not be requested, thus,', &
       ' rci_request will not take the value 2'
      endif
      write( *,'(a)') '+++'
      if (ipar(11).ne.0) then
      	write(*,'(a,i1,a,a)') 'as ipar(11)=',ipar(11),', the', &
       ' preconditioned fgmres iterations will be performed, thus,'
      	write(*,'(a,a)') 'the preconditioner action will be requested', &
       ' via rci_request=3'
      else
      	write(*,'(a,i1,a,a)') 'as ipar(11)=',ipar(11),', the', &
       ' preconditioned fgmres iterations will not be performed,'
      	write( *,'(a)') 'thus, rci_request will not take the value 3'
      endif
      write( *,'(a)') '+++'
      if (ipar(12).ne.0) then
      	write(*,'(a,i1,a,a)')'as ipar(12)=',ipar(12),', the automatic', &
       ' test for the norm of the next generated vector is not'
      	write( *,'(a,a)') ' equal to zero up to rounding and', &
       ' computational errors will be performed,'
      	write( *,'(a)') 'thus, rci_request will not take the value 4'
      else
      	write(*,'(a,i1,a,a)')'as ipar(12)=',ipar(12),', the automatic', &
       ' test for the norm of the next generated vector is'
      	write(*,'(a,a)') 'not equal to zero up to rounding and', &
      	' computational errors will be skipped,'
      	write(*,'(a,a)') 'thus, the user-defined test will be requested', &
       ' via rci_request=4'
      endif
      write( *,'(a)') '+++'

end subroutine mkl_gmres_print_info

!==========================================================================================

end module iter_solvers
