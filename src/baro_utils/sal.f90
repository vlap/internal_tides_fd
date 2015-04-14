module sal

     use precisions, only: wp, cwp
     use my_trigs
     use my_sparse
     use my_sparse_aggregate
     use control
     use generate_grid
     use save_load
     use interpolate
     use spherepack_iface
	 use dispmodule

contains

!**********************************************************************************************************
!**********************************************************************************************************
subroutine calc_sal(ncpts, cpts, P, GD, dir_grid, dir_cols, dir_sols, beta0)

!% [hsal_col,beta0,beta_col]=calc_hsal_col_merc(h_col,flags)
!% Given a column vector of surface heights h_col, distributed
!% globally according to information in , calculates the
!% corresponding h_SAL in column format, along with a local
!% beta term (h_SAL approx = beta.*h), and a global average beta0.

implicit none

     character(len=*),intent(in)    :: cpts
     integer, intent(in)			:: ncpts
     character(len=2)  				:: cpt
     type(params), intent(in)		:: P
     type(grid_dims), intent(in)	:: GD

	 complex(cwp),allocatable	 :: h(:)
     integer	    		::  nh, ccpt

     real(wp)			   :: beta0(ncpts)
     real(wp), allocatable :: beta(:)
     complex(cwp) :: beta0_cmplx
     complex(cwp), allocatable :: hsal(:), beta_cmplx(:)

     integer		:: j, istat
!     real    ::      T1, T2
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     character(len = *) :: dir_grid, dir_cols, dir_sols

!***************************************************

	if (P%messages > 0) then
		print *, ""
		write(*, '("Calculate hsal and beta: ")')!, advance = 'no')
	endif

	do ccpt = 1, ncpts

	    cpt=cpts(2*ccpt-1:2*ccpt)

		if (P%messages >= 1) then
			write(*, '(i2, ". for ", a, ": ")', advance='no') ccpt, cpt
		endif

		call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat', nh)

	    if (P%sal_scheme == 2) then
			call calc_hsal_col(h, P, GD, dir_cols, dir_grid, hsal, beta0_cmplx)
	      	beta0(ccpt) = real(beta0_cmplx, wp)
	     	call save_vector(hsal, dir_sols //'temp/'// cpt // '_hsal' // '.dat')
	     	deallocate(hsal)
	    elseif (P%sal_scheme == 3) then
	      	call calc_hsal_col(h, P, GD, dir_cols, dir_grid, hsal, beta0_cmplx, beta_cmplx)
	      	beta0(ccpt) = real(beta0_cmplx, wp)
	      	allocate(beta(nh), stat = istat)
	      	beta = real(beta_cmplx, wp)
	      	where (beta<P%betamin) beta = P%betamin
	      	where (beta>P%betamax) beta = P%betamax
	     	call save_vector(hsal, dir_sols //'temp/'// cpt // '_hsal' // '.dat')
	     	call save_vector(beta, dir_sols //'temp/'// cpt // '_betasal' // '.dat')
	     	deallocate(hsal, beta_cmplx, beta)
	    endif

    	deallocate(h)
	enddo

! Clean up: deallocate work arrays used in SHT routines
	call cleanup('REG')

end subroutine calc_sal

!**********************************************************************************************************
!**********************************************************************************************************
subroutine calc_hsal_col(h_col, P, GD, dir_cols, dir_grid, hsal_col, beta0, beta_col)

!% [hsal_col,beta0,beta_col]=calc_hsal_col_merc(h_col,flags)
!% Given a column vector of surface heights h_col, distributed
!% globally according to information in , calculates the
!% corresponding h_SAL in column format, along with a local
!% beta term (h_SAL approx = beta.*h), and a global average beta0.

implicit none

     complex(cwp),intent(in):: h_col(:)
     integer, allocatable   :: hp(:,:)
     integer	    :: nph, nth, nh
     type(params)	:: P
     type(grid_dims):: GD

     complex(cwp), allocatable, optional:: beta_col(:)
     complex(cwp), optional				:: beta0

     type (triplet_cmplx)	      :: triplet_h
     complex(cwp),allocatable :: h(:,:)

!     real(wp), pointer    :: th_h(:), ph_h(:)
     complex(cwp), allocatable	:: hsal(:,:), beta(:,:), hsal_col(:)
     integer		:: j, istat
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     character(len = *) :: dir_grid, dir_cols

!***************************************************

!***************************************************
!	PERFORMANCE: START
!***************************************************
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )


! shortcuts
nph = GD%nph
nth = GD%nta
nh = GD%nh
! load hp and project hcol to the grid
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

    call init_sparse_vals(triplet_h, hp(:,1),hp(:,2),h_col, nph,nth, nh)
    call coo2full(triplet_h, h)
    call dealloc_sparse(triplet_h)

!	call save_matrix(h, dir_grid // 'm2' // '_h_grid' // '.dat')

	if (present(beta_col)) then ! get hsal, beta0, beta
		call calc_hsal_gridded(h, dir_grid, P, GD, hsal, beta0, beta)
    elseif (present(beta0)) then ! get hsal, beta0
		call calc_hsal_gridded(h, dir_grid, P, GD, hsal, beta0)
	else ! get hsal
		call calc_hsal_gridded(h, dir_grid, P, GD, hsal)
    end if

!print *, count(isnan(abs(beta)))

	where (abs(h)==0) hsal = 0
!	call save_matrix(hsal, dir_grid // 'm2' // '_hsal_grid' // '.dat')
!	call save_matrix(beta, dir_grid // 'm2' // '_beta_grid' // '.dat')

 ! write hsal and beta in columns
    allocate(hsal_col(nh), stat = istat)
    hsal_col = 0
    do j = 1, nh
		hsal_col(j) = hsal(hp(j,1),hp(j,2))
	end do

	if (present(beta_col)) then
	    allocate(beta_col(nh), stat = istat)
	    do j = 1, nh
			beta_col(j) = beta(hp(j,1),hp(j,2))
		end do
    end if

!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )
if (P%messages >= 1) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
!	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall:'&
!							//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')

elseif (P%messages == 0) then
	write(*, '("")')
endif

end subroutine calc_hsal_col

!**********************************************************************************************************
!**********************************************************************************************************

subroutine calc_hsal_gridded(h, dir_grid, P, GD, hsal, beta0, beta)

!% Given h on a (phi,tau) grid, calculates hsal and beta
!% on the corresponding grid. Also calculates a single
!% globally integrated value beta0.

implicit none

!     integer, allocatable   :: hp(:,:)
     complex(cwp), allocatable	:: h(:,:), h_sht(:,:)
     complex(cwp), allocatable	:: hsal(:,:), hsal_sht(:,:)
     complex(cwp), allocatable	:: tmp1(:,:), tmp2(:,:)!, beta_tmp(:,:)
     real(wp), allocatable	:: tmp3(:,:), tmp4(:,:)

     complex(cwp), allocatable, optional:: beta(:,:)
     complex(cwp), optional				:: beta0

     type(params)		:: P
     type(grid_dims)	:: GD

     integer			  :: nph, nth, nph_sht, nth_sht ! , nh
     real(wp)             :: rhoe, rhoo, dph_sht, dth_sht

     real(wp), allocatable	:: th_hg(:), ta_ug(:), ph_hg(:), ph_sht(:), th_sht(:), metrics(:, :)

!     real(wp), allocatable	:: test_r(:,:), test_i(:,:)
!     complex(cwp), allocatable	:: test_spec(:), test(:,:), test_int(:,:)

     integer, allocatable :: mask(:,:)!, frontmat_ph_f(:,:), backmat_th_f(:,:)
     real(wp), allocatable :: H_hg(:, :)
     integer            :: j, c, istat, ntrunc, nf, nnz, ncurrent

     type (triplet_int)	:: frontmat_ph, backmat_th

     character(len = *) :: dir_grid!, dir_cols

     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max

!==========================================================================================
! DISP MODULE SETTINGS
!==========================================================================================
call tostring_set(rfmt='F12.1')

! load ph_vg, th_ug
     call load_alloc_vector(ph_hg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(th_hg, dir_grid // 'th_ug.dat')

! shortcuts
nph = GD%nph
nth = GD%nta
!nh =  GD%nh

rhoe = P%rhoe
rhoo = P%rhoo

!***************************************************
! Define truncation cut-off
! On a regular grid spanning -pi/2<th<pi/2 and 0<ph<2*pi
! the condition (ntrunc <= nth-1) is equiv. to (2*ntrunc <= nph)

if (nph > 2*P%ntrunc) then
	ntrunc = P%ntrunc
else
	ntrunc = nph/2
endif
!print *, "ntrunc", ntrunc
!***************************************************
! Construct a regular spaced lon-lat grid for the SHT routine
! as the basis take the ph_hg partition of [0, 2*pi) interval (ph_sht = ph_hg)
! Lat coord th_sht must span [-pi/2, pi/2] and must start at the north pole

!***************************************************
!	PERFORMANCE: START
!***************************************************
if (P%messages >= 3) then
	write(*, *) ""
	write(*, '("   a) Interpolate h to the reg ")', advance = 'no')
	call disp('['//tostring(2*P%ntrunc)// ' x ' // tostring(P%ntrunc) //'] lat-lon grid...', advance = 'no')
endif

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )
!***************************************************

!nph_sht = nph
nph_sht = 2*ntrunc
dph_sht = 2*pi/nph_sht

nth_sht = nph_sht/2 + 1
dth_sht = pi/(nth_sht-1)

if (dph_sht /= dth_sht) then
	print *, 'dph_sht /= dth_sht in SHT subroutine calc_hsal_gridded' ! implies that nph is even!!!!!!!!!!!!
	stop
end if

! initialises the sht grid structure
    allocate(ph_sht(nph_sht), stat = istat)
    ph_sht = (/ ( (2*pi/nph_sht*j), j=0,nph_sht-1 ) /)	! the grid for the SHT routine spans longitude points
			 											! phi(j) = (j-1)*2*pi/nlon
! Th must be the colatitude needed by the SHT routine (increases from the North pole southwards)
    allocate(th_sht(nth_sht), stat = istat)
    th_sht = (/ ( (pi/2 - pi/(nth_sht-1)*j), j=0,nth_sht-1 ) /)


! interpolate h to the grid
    call bilinear_2d(nph, nth, ph_hg, th_hg, h, nph_sht, nth_sht, ph_sht, th_sht, h_sht, P%lib)!, dir_grid)
!      call save_topo(dir_grid // 'h_orig.dat', nph, nth, ph_hg, th_hg, h)
!      call save_topo(dir_grid // 'h_interp.dat', nph_sht, nth_sht, ph_sht, th_sht, abs(h_sht))

!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (P%messages >= 3) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
!	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall:'&
!							//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')
elseif (P%messages == 0) then
	write(*, '("")')
endif
!***************************************************

! do the main calculation to find hsal on the grid
	if (P%save_sht==0) then
		call calc_hsal('0', nph_sht, nth_sht, h_sht, rhoe, rhoo, ntrunc, hsal_sht, (P%messages >= 3))
	else
		call calc_hsal(dir_grid // 'temp/', nph_sht, nth_sht, h_sht, rhoe, rhoo, ntrunc, hsal_sht, (P%messages >= 3))
	endif
!  	call save_topo(dir_grid // 'h_sal.dat', nph_sht, nth_sht, ph_sht, th_sht, abs(hsal_sht))

!print *, 'beta wo mask', sum( hsal_sht*conjg(h_sht) ) / sum( h_sht*conjg(h_sht) )
!print *, 'beta abs', sum( abs(hsal_sht)*abs(h_sht) ) / sum( abs(h_sht)**2 )
!    deallocate(h_sht)
!***************************************************
!	PERFORMANCE: START
!***************************************************
if (P%messages >= 3) then
	write(*, '("   d) Interpolate hsal onto the original grid... ")', advance = 'no')
endif

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )
!***************************************************
! interpolate hsal back to the original grid
    call bilinear_2d(nph_sht, nth_sht, ph_sht, th_sht, hsal_sht, nph, nth, ph_hg, th_hg, hsal, P%lib)!, dir_grid)
!      	call save_topo(dir_grid // 'h_sal_interp.dat', nph, nth, ph_hg, th_hg, hsal)
!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (P%messages >= 3) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
!	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall:'&
!							//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')
elseif (P%messages == 0) then
	write(*, '("")')
endif
!***************************************************

!    deallocate(ph_sht, th_sht)
    deallocate(hsal_sht)

!==========================================================================================
!***************************************************
! TEST SH TRANSFORM
!print *, "Test SHT start"
! comment:     deallocate(ph_sht, th_sht)
!
!!allocate(test(nph_sht, nth_sht),test_r(nph_sht, nth_sht),test_i(nph_sht, nth_sht), stat = istat)
!!allocate(test_spec((ntrunc+1)*(ntrunc+2)/2), stat = istat)
!!	do j=1, nth_sht
!!		test(:,j) = exp(i*ph_sht(:))*cos(th_sht(j))
!!	end do
!
!    call save_topo(dir_grid // 'h_test.dat', nph_sht, nth_sht, ph_sht, th_sht, h_sht)
!	test_r = real(h_sht,wp)
!	test_i = imag(h_sht)
!
!	call grdtospec(test_r,test_spec,'0', 'REG')
!	call spectogrd(test_spec, test_r,'0', 'REG')
!
!	call grdtospec(test_i,test_spec,'0', 'REG')
!	call spectogrd(test_spec, test_i, '0', 'REG')
!
!	test = cmplx(test_r, test_i, kind = cwp)
!	call save_topo(dir_grid // 'h_test_inv.dat', nph_sht, nth_sht, ph_sht, th_sht, test)
!
!print *, "Test SHT end"
!stop
!***************************************************

!***************************************************
! TEST bilinear interp
!print *, "Test bilinear start"
!
!allocate(test(nph, nth),test_int(nph_sht, nth_sht), stat = istat)
!	do j=1, nth
!		test(:,j) = exp(i*ph_hg(:))*cos(th_hg(j))
!	end do
!    call save_topo(dir_grid // 'h_test.dat', nph, nth,  ph_hg, th_hg, test)
!
!	call bilinear_2d(nph, nth, ph_hg, th_hg, test, nph_sht, nth_sht, ph_sht, th_sht, test_int, P%lib)
!	call save_topo(dir_grid // 'test_int.dat', nph_sht, nth_sht, ph_sht, th_sht, test_int)
!
!	call bilinear_2d(nph_sht, nth_sht, ph_sht, th_sht, test_int, nph, nth, ph_hg, th_hg, test, P%lib)
!	call save_topo(dir_grid // 'h_test_back.dat', nph, nth,  ph_hg, th_hg, test)
!
!print *, "Test bilinear end"
!***************************************************
!==========================================================================================


! Now calculate beta0, beta for the iterative routine
    if (present(beta0)) then

!***************************************************
!	PERFORMANCE: START
!***************************************************
if (P%messages >= 3) then
	write(*, '(a)', advance='no') '   e) Calculate beta (avrgd over '
	call disp(tostring(P%sal_avg)//'º)...', advance='no')
endif

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )
!***************************************************

	     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
	     call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')

	    ! Mask out the non-ocean hsal
	    allocate(mask(nph, nth), stat = istat)
	    mask = 0
	    where (H_hg > 0) mask = 1

	    allocate(metrics(nph, nth), stat = istat)
		    if (P%coor == 1) then
	    	do j = 1,nth
	    		metrics(:, j) = 1 / cosh(ta_ug(j))**2 ! Include metrics for Mercator coords
	    	end do
	    elseif (P%coor == 2) then
	    	do j = 1,nth
	    		metrics(:, j) = cos(th_hg(j)) ! Include metrics for Lat-Lon coords
	    	end do
	    endif

	    deallocate(ta_ug, H_hg)

!	call save_topo(dir_grid // 'mask.dat', nph, nth, ph_hg, th_hg, mask)
!	call save_matrix(metrics, dir_grid // 'tmp' // '.dat')
		beta0 = sum( mask*metrics*hsal*conjg(h) ) / sum( mask*metrics*h*conjg(h) )

!print *, 'beta0', beta0

	if (present(beta)) then

		allocate(beta(nph, nth), stat = istat)
		beta = beta0

		if (nph*P%sal_avg/360 > 1.5) then
!			averaging in ph
			nf=ceiling( nph*P%sal_avg/360 ) ! will average over 2*nf+1 cells (extra nf in ph and th directions)
			nnz=nph*(2*nf+1)
    		call init_sparse_0(frontmat_ph, nph, nph, nnz)

			ncurrent = 1
			do c = 1,nph
				frontmat_ph%indi(ncurrent : ncurrent + 2*nf) = c
				frontmat_ph%indj(ncurrent : ncurrent + 2*nf) = 1 + modulo( (/(j, j=c-nf,c+nf)/)-1, nph)
				frontmat_ph%vals(ncurrent : ncurrent + 2*nf) = 1

				ncurrent = ncurrent + 2*nf+1
			end do

!			averaging in th
!			first count the number of elements
			nnz = 0
			do c = 1,nth
				if (P%coor == 1) then
					nf=ceiling( nph*P%sal_avg/360 /cos(th_hg(c)) ) ! 1/cos(th_ug(cta)) bc of Mercator coords
				else
					nf=ceiling( nph*P%sal_avg/360 )
				endif
				nnz = nnz + min(c+nf,nth) - max(1,c-nf) + 1
			end do
			call init_sparse_0(backmat_th, nth, nth, nnz)

			ncurrent = 1
			do c = 1,nth
				if (P%coor == 1) then
					nf=ceiling( nph*P%sal_avg/360 /cos(th_hg(c)) ) ! 1/cos(th_ug(cta)) bc of Mercator coords
				else
					nf=ceiling( nph*P%sal_avg/360 )
				endif
				nnz = min(c+nf,nth) - max(1,c-nf) + 1

				backmat_th%indi(ncurrent : ncurrent + (nnz-1)) = (/(j, j=max(1,c-nf),min(c+nf,nth))/)
				backmat_th%indj(ncurrent : ncurrent + (nnz-1)) = c
				backmat_th%vals(ncurrent : ncurrent + (nnz-1)) = 1

				ncurrent = ncurrent + nnz
			end do

!			allocate(beta(nph, nth), stat = istat)
!***************************************************
!      	call save_topo(dir_grid // 'tmp.dat', nph, nth, ph_hg, th_hg, abs(hsal*conjg(h)))
!      	call save_topo(dir_grid // 'tmp_1.dat', nph, nth, ph_hg, th_hg, abs(matmul( hsal*conjg(h),backmat_th )))
!      	call save_topo(dir_grid // 'tmp_2.dat', nph, nth, ph_hg, th_hg, &
!      					abs(matmul(matmul(frontmat_ph, hsal*conjg(h) ), backmat_th)))


!	Calculate using multiplication of full matrices (slow x2)

!    call coo2full(frontmat_ph, frontmat_ph_f)
!    call coo2full(backmat_th, backmat_th_f)
!			!      call save_sparse(frontmat_ph, dir_grid // 'frontmat_ph.dat')
!			!      call save_sparse(backmat_th, dir_grid // 'backmat_th.dat')
!			!      call save_matrix(frontmat_ph_f, dir_grid // 'frontmat_ph_f.dat')
!			!      call save_matrix(backmat_th_f, dir_grid // 'backmat_th_f.dat')
!  call CPU_Time(T1)
!***************************************************
!		beta = mask*matmul(matmul(frontmat_ph_f, metrics*hsal*conjg(h) ), backmat_th_f) / &
!			        matmul(matmul(frontmat_ph_f, metrics*h*conjg(h) ),    backmat_th_f)
!***************************************************
!      	call save_topo(dir_grid // 'beta.dat', nph, nth, ph_hg, th_hg, abs(beta))
!  call CPU_Time(T2)
!  call disp ('Done in ' // tostring(T2-T1) // 's CPU')

!	CHECK PERFORMANCE
!  call CPU_Time(T1)
			allocate(tmp1(nph, nth), tmp2(nph, nth), stat = istat) ! complex valued
			call coo_mat_mul(frontmat_ph, nph, nth, metrics*hsal*conjg(h), P%lib, tmp1)
!		   	call save_topo(dir_grid // 'tmp1.dat', nph, nth, ph_hg, th_hg, abs(tmp1))
			call mat_coo_mul(nph, nth, tmp1,backmat_th, P%lib, tmp2)
!      		call save_topo(dir_grid // 'tmp2.dat', nph, nth, ph_hg, th_hg, abs(tmp2))
			beta = tmp2
			deallocate(tmp1, tmp2)

			allocate(tmp3(nph, nth),tmp4(nph, nth), stat = istat) ! real valued
			call coo_mat_mul(frontmat_ph, nph, nth, metrics*(abs(h)**2), P%lib, tmp3)
			call mat_coo_mul(nph, nth, tmp3,backmat_th, P%lib, tmp4)
			where (mask .ne. 0)
				beta = beta / tmp4
			elsewhere
				beta = 0
			end where
			deallocate(tmp3, tmp4)
!  call CPU_Time(T2)
!	        call disp ('Done in ' // tostring(T2-T1) // 's CPU')
!      		call save_topo(dir_grid // 'beta2.dat', nph, nth, ph_hg, th_hg, abs(beta))
		else
			where (mask .ne. 0)
				beta = hsal / h ! beta = mask*metrics*hsal*conjg(h) / mask*metrics*h*conjg(h)
			elsewhere
				beta = 0
			end where
		endif

	end if
!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (P%messages >= 3) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
elseif (P%messages == 0) then
	write(*, '("")')
endif
!***************************************************

  end if


call tostring_set_factory()

!    deallocate(h)

end subroutine calc_hsal_gridded
!**********************************************************************************************************

subroutine calc_hsal(save_sht, nph, nth, h, rhoe, rhoo, ntrunc, hsal, messages)! ph, th,

!% returns hsal globally, so might want to
!% pre-multiply by the ocean-mask...
!% h must be given on the grid defined by shtmats

implicit none

     character(len = *), parameter :: gridtype = 'REG' ! (regularly spaced grid)
     character(len = *), intent(in) :: save_sht

     integer, intent(in)		:: nph, nth, ntrunc
     integer					:: nmdim ! dimension((ntrunc+1)*(ntrunc+2)/2) where mtrunc = nth - 1
     complex(cwp), intent(in)	:: h(nph, nth)
     real(wp), intent(in)		:: rhoe, rhoo

     real(wp)		:: hsal_r(nph, nth), hsal_i(nph, nth)
     complex(cwp), allocatable	:: hsal(:, :)

     complex(cwp), dimension((ntrunc+1)*(ntrunc+2)/2)	:: hspec_r, hspec_i, hsalspec_r, hsalspec_i
     real(wp), dimension((ntrunc+1)*(ntrunc+2)/2)     :: myf

!     real(wp)			:: th(:), ph(:)
     integer            :: j, m,n, istat, nvals(nth)
	 real(wp)			:: hn(nth), kn(nth)
	 logical			:: f_exist

	 !	Some useful indices for addressing SHT coefficients
	 integer, dimension(nth*(nth+1)/2) :: indxn!, indxm

	 logical :: messages
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max

!	some indexing short
	nmdim = (ntrunc+1)*(ntrunc+2)/2! dimension((ntrunc+1)*(ntrunc+2)/2) where mtrunc = nth - 1
!	indxm = (/((m,n=m,ntrunc),m=0,ntrunc)/)
	indxn = (/((n,n=m,ntrunc),m=0,ntrunc)/)


!***************************************************
!	PERFORMANCE: START
!***************************************************
if (messages) then
	write(*, '(a)', advance='no') '   b) Spherical harmonic analysis... '
endif

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )
!***************************************************
! Apply SHT to convert gridded input array (datagrid) to complex spectral coefficients
!	Do separately for Im and Re parts
		call grdtospec(real(h,wp),hspec_r,save_sht, gridtype)
		call grdtospec(imag(h),hspec_i, save_sht,gridtype)

!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (messages) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
endif
!***************************************************

!	Get love numbers
!print *, "love numbers"
inquire( file=save_sht//'love_num_h.dat', exist=f_exist )
if ( (save_sht .eq. '0') .or. (.not. f_exist ) ) then
		nvals  = (/ ( (j), j=1,nth ) /)
		call get_load_love(nth, nvals, hn, kn)
	if (save_sht .ne. '0') then
		call save_vector(hn, save_sht // 'love_num_h.dat')
		call save_vector(kn, save_sht // 'love_num_k.dat')
	endif
else
		call load_vector(hn, save_sht // 'love_num_h.dat')
		call load_vector(kn, save_sht // 'love_num_k.dat')
endif
!	Calculate the SAL coefficients
	myf(1) = 0
	myf(2:nmdim) = (1 + kn(indxn(2:nmdim)) - hn(indxn(2:nmdim))) / (2*indxn(2:nmdim) + 1)

!	Calculate the SAL in spherical harmonics
	hsalspec_r = (3*rhoo/rhoe)*myf*hspec_r
	hsalspec_i = (3*rhoo/rhoe)*myf*hspec_i

!***************************************************
!	PERFORMANCE: START
!***************************************************
if (messages) then
	write(*, '(a)', advance='no') '   c) Spherical harmonic synthesis... '
endif

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )
!***************************************************
! 	Inverse transform to obtain hsal on the grid
!    allocate(hsal(nph, nth), stat = istat)
!	call spectogrd(hsalspec,hsal,gridtype)
		call spectogrd(hsalspec_r,hsal_r, save_sht,gridtype)
		call spectogrd(hsalspec_i,hsal_i, save_sht,gridtype)
!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (messages) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
endif
!***************************************************
!	call save_matrix(hsal_r, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/' //&
!								 'm2' // '_hsal_r_grid' // '.dat')
!	call save_matrix(hsal_i, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/' //&
!								 'm2' // '_hsal_i_grid' // '.dat')

!	Export resulting hsal
	allocate(hsal(nph, nth), stat = istat)
	hsal = cmplx(hsal_r, hsal_i, kind=cwp)

end subroutine calc_hsal

!**********************************************************************************************************
!**********************************************************************************************************

subroutine calc_heq_h(heq, cpt, nh, latP, lonP, coor, dir_cols, dir_grid)

!% calculates the (complex) equilibrium tide heq
!% on the h grid, for the specified domain
implicit none

     real(wp), allocatable        :: ta_h(:), ph_vg(:)!, ta_ug(:)
     integer, allocatable       :: hp(:,:)

     real(wp), intent(in)     :: latP, lonP
     integer, intent(in)    :: nh, coor

     character(len=2)   :: cpt
     real(wp)             :: th0, ph0

     real(wp), pointer    :: th_h(:), ph_h(:)
     integer            :: status1, status2, status3!, j

     character(len = *) :: dir_grid, dir_cols

     complex(cwp), allocatable, dimension(:) :: heq

! load hp, ph_vg, ta_h
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
!     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')

!%===============
!% calculate heq
!%===============
    if (associated(th_h)) deallocate(th_h)
    allocate(th_h(nh), stat = status1)
    if (associated(ph_h)) deallocate(ph_h)
    allocate(ph_h(nh), stat = status2)

    allocate(heq(nh), stat = status3)

!th_h = (/ ( (ta2th(ta_h(j))), j=1,nh ) /)
th_h = ta2th(ta_h, coor)
!th_h = (/ ( (ta2th(ta_ug(hp(j,2)))), j=1,nh ) /)
ph_h = ph_vg(hp(:,1)) ! 1 for ph coordinate

     deallocate(ta_h, ph_vg, hp)

th0 = d2r(latP)
ph0 = d2r(lonP)

 call calc_heq_i(heq, cpt, th_h, ph_h, th0,ph0)

!allocate(h_abs(size(ph_vg),size(ta_ug)), stat = status3)
!h_abs = 0
!
!do j = 1,size(th_h)
!    h_abs(hp(j,1), hp(j,2)) = abs(heq(j))
!enddo

!call save_vector(ta_ug, '/home/amsta/amtvl/Desktop/scratch/Tides/data/cols_ta.dat')
!call save_vector(ph_vg, '/home/amsta/amtvl/Desktop/scratch/Tides/data/cols_ph.dat')
!call save_matrix(h_abs, '/home/amsta/amtvl/Desktop/scratch/Tides/data/mats_h_abs.dat')

end subroutine calc_heq_h

!**********************************************************************************************************
!**********************************************************************************************************

subroutine calc_heq_i(heq, cpt,th,ph,th0,ph0)
!
!% returns a complex equilibrium tide heq, where the actual
!% equilibrium tide is Re ( heq(th,ph) e(-i w t) ).
!% if th and ph are the same size, returns heq at that size.
!% if th is column and ph is row, returns heq as th x ph matrices.
!% th0 and ph0 are optional values by which the grid has been
!% rotated from geo coordinates.
implicit none

     character(len=2)   :: cpt
     type (tide_params) :: pars
     real(wp)           :: th0, ph0
     real(wp), pointer  :: th(:), ph(:)

     complex(cwp), allocatable, dimension(:) :: heq(:)

if (size(th)/=size(ph)) then
    write(*, '("** See function calc_heq_i")')
    stop
    return
endif

     pars = get_pars(cpt)

if ( (cpt == 'k2') .or. (cpt == 'm2') .or. (cpt == 'n2') .or. (cpt == 's2') .or. (cpt == 't2') ) then

    heq = pars%amp*pars%lovef*exp(-2*i*ph0) &
                *( cos(th0)*cos(th)*cos(ph)+sin(th0)*sin(th)-i*cos(th)*sin(ph) )**2

elseif ( (cpt == 'k1') .or. (cpt == 'p1') .or. (cpt == 'q1') .or. (cpt == 'o1') .or. (cpt == 't1') ) then

    heq = pars%amp*pars%lovef*exp(-i*ph0) &
                *( cos(2*th0)*sin(2*th)*cos(ph) + sin(2*th0)*(sin(th)**2) &
                  -sin(2*th0)*(cos(th)**2)*(cos(ph)**2) + i*sin(th0)*(cos(th)**2)*sin(2*ph) &
                  -i*cos(th0)*sin(2*th)*sin(ph) )

elseif ( (cpt == 'mm') .or. (cpt == 'mf') ) then

    heq = pars%amp*pars%lovef*(0.5 - 1.5*(cos(th0)*sin(th)-sin(th0)*cos(th)*cos(ph))**2)

endif


end subroutine calc_heq_i

!**********************************************************************************************************

end module sal
