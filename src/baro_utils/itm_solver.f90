module itm_solver_mod

!     use my_trigs
!     use generate_global_matrices
!	 use generate_grid
	 use control
	 use interpolate
     use unsym_solvers
     use iter_solvers
     use my_sparse_aggregate
     use my_sparse
     use sal
     use itd
     use save_load
     use dispmodule
     use precisions, only: wp, cwp

     contains

!==========================================================================================
subroutine forced_itm_solver(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

    implicit none
!%   Calculates internal tide field forced by given barotropic tide (no iterations).
     character(len=*) :: cpts
     type(params)	:: P
     type(grid_dims):: GD
     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols
!==========================================================================================


		call itm_solver_linear(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)
		call itm_drag_m0(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)
	!==========================================================================================
	!      %==============================
	!      % plot global tides and errors
	!      %==============================

!		if  ( P%graphics == 2 ) then
!	!		Use MATLAB functions to export the data in M-files into binaries
!			call system('xterm -e matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
!			'"addpath(genpath(''' // matlab_dir // ''')); baro_v1_plot_conv(''' // &
!			trim(save_gridid)//''', '// tostring(itn) // ', ''' // cpts // ''', '//tostring(P%load_sol)//'); waitforbuttonpress; exit;" ')
!		end if
	!==========================================================================================

!=========================================================================================
!=================================DONE====================================================
!=========================================================================================
write(*, '("========================================================================")')
write(*, '("========================================================================")')

end subroutine forced_itm_solver
!==========================================================================================
subroutine itm_solver_linear(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

    implicit none
!%   Calculates a baroclinic tides. Uses u/v/h + volume transport formulation
!%
     integer :: nph, nth
     integer :: nu, nv, nh, np
!     real(wp):: latP, lonP
     real(wp):: g, re!, cdg, Q, beta
     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols
     type (tide_params) :: pars
     type(params)	:: P
     type(grid_dims):: GD

	 logical :: 	 CGS_ready = .false.
     integer ::      j, istat, iter
     real    ::      T1_itn, T2_itn, T1_cpt, T2_cpt ! for measuring CPU (NOT REAL TIME!)
     integer :: 	 wall_t1_itn, wall_t2_itn, wall_t1_cpt, wall_t2_cpt, clock_rate, clock_max
     integer, allocatable :: hp(:, :),up(:, :),vp(:, :)

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt, nmodes, cmode
     integer :: iu(GD%nu), iv(GD%nv), ih(GD%nh)

     real(wp), allocatable	:: issponge_u(:, :), issponge_v(:, :), issponge_h(:, :)
     complex(cwp), allocatable :: temp_vec(:)

     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr) :: mat_csr
     type (csr_cmplx) :: mat_csr_cmplx
     type (csc_cmplx) :: mat_csc_cmplx

!	 type (csr) :: tmp_csr
!	 type (csr_cmplx) :: tmp_csr_cmplx

!	 type (triplet) :: tempmat
!     integer ::      nmat, filled

!	 type (triplet) :: DUU, DVU, DUV, DVV
!     complex(cwp), allocatable :: Qu(:), Qv(:)!, tmp(:), Du(:), Dv(:)
!	 real(wp), allocatable :: KU(:), KV(:)
!     integer, allocatable, dimension(:) :: ufactor,vfactor

     type (csr_cmplx)			:: speye
     integer, allocatable		:: bcdiag(:)
     complex(cwp), allocatable	:: rhs(:,:), rhs_cur(:)

     complex(cwp), allocatable :: uvp(:), u(:), v(:), pr(:)

!==========================================================================================
! SHORTCUTS
!==========================================================================================
!latP = P%latP
!lonP = P%lonP
g = P%g
re = P%re
nmodes = P%n_modes

nph = GD%nph
nth = GD%nta

nu = GD%nu
nv = GD%nv
nh = GD%nh
np = GD%np

!==========================================================================================
write(*, '("------------ INTERNAL TIDES ------------")')


ncpts=len(cpts)/2

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (nu + j), j=1,nv ) /)
     ih = (/ ( (nu + nv + j), j=1,nh ) /)

!%==========================
!% solve for each component
!%==========================
do ccpt = 1, ncpts

    cpt=cpts(2*ccpt-1:2*ccpt)
    pars = get_pars(cpt)

		write(*, '(i2, ") Calc barotropic forcing for ", a, "... ")', advance='no') ccpt, cpt
		!%=======================
		!% calculate rhs forcing
		!%=======================
		! calc rhs for all nmodes
		call itm_rhs(rhs, cpt, P, GD, nmodes, dir_cols, dir_grid, dir_mats, dir_sols); write(*, '(": ")')

!		call save_vector(rhs, dir_cols //'temp/'// 'itm_rhs_0.dat')
		!%=============
		!% impose bcs:
		!%=============
		call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
		do cmode = 1, nmodes
			rhs(:,cmode) = rhs(:,cmode) * bcdiag
		enddo
		!call save_vector(bcdiag, dir_cols // 'bcdiag' //'.dat')
		deallocate(bcdiag)

		! load sponge
		call load_alloc_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//cpt(2:2)//'.dat')
		call load_alloc_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//cpt(2:2)//'.dat')
		call load_alloc_matrix(issponge_h, dir_grid//'itm/'//'issponge_h'//cpt(2:2)//'.dat')

	do cmode = 1, nmodes

		write(*, '("mode ",i2, ": ")', advance='no') cmode
		!%==================================================
		!% Load the basic sparse matrix generated by write_itm_mats
		!% (includes all non-diag terms: rotation, c.o.mass, grad(P))
		!%==================================================
		call load_alloc_sparse(mat_csr, dir_mats // 'itm/temp/' // 'itm_mat_csr_init_'//&
	     								cpt(2:2)//'_'//tostring(cmode)//'.dat')
		!  now transform into cmplx matrix
		call init_sparse_vals(mat_csr_cmplx,mat_csr%indi,mat_csr%indj,mat_csr%vals, &
		                                mat_csr%ni, mat_csr%nj, mat_csr%nz)
		call dealloc_sparse(mat_csr)
		!  %=============
		!  % impose bcs:
		!  %=============
		call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
		!mat=bcmat*mat; !fill the row of mat with zeroes for boundary variables (then will get -i*omega)
		!rhs=bcmat*rhs;
		call dia_csr_mul(real(bcdiag,wp), mat_csr_cmplx, skit)
		deallocate(bcdiag)

		!%===================================
		!	Diagonal -i*omega AND sponge damping terms
		!%===================================
		allocate(temp_vec(np), stat=istat)
			temp_vec(iu) = pars%omega0*(cmplx(P%bulk_damp, -1, kind=cwp) + P%sponge_damp*issponge_u(:, cmode))
			temp_vec(iv) = pars%omega0*(cmplx(P%bulk_damp, -1, kind=cwp) + P%sponge_damp*issponge_v(:, cmode))
!			temp_vec(ih) = pars%omega0*(cmplx(P%bulk_damp, -1, kind=cwp) + P%sponge_damp*issponge_h(:, cmode))
			temp_vec(ih) = pars%omega0*cmplx(0, -1, kind=cwp)

		!     mat = mat + (-i*omega + omega*P%sponge_damp*issponge)*speye(np)
	      speye = init_speye(np, temp_vec)
	      call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), speye, mkl)
	      deallocate(temp_vec)
	      call dealloc_sparse(speye)

	!%=========================================================================================
	!write(*, '("DONE")')
	!==========================================================================================
		!     Save the final version of the matrix
		if (P%messages>=3) then
		      call csr2coo(mat_csr_cmplx, mat_coo_cmplx, P%lib)
		      call save_sparse(mat_coo_cmplx, dir_mats // 'mat_' // cpt //'.dat')
		      call dealloc_sparse(mat_coo_cmplx)
		endif

	!  %=======
	!  % solve
	!  %=======
		if (P%messages >= 1) then
			call CPU_Time(T1_cpt)
			call system_clock ( wall_t1_cpt, clock_rate, clock_max )
		endif

	!********************************************************************
		if (P%solver .eq. 'umfpack') then
		!	USE UMFPACK TO SOLVE THE SYSTEM

		      ! convert from csr to csc for umfpack solver (csc is convert to 0-based array)
		        call csr2csc(mat_csr_cmplx, mat_csc_cmplx, P%lib)
		        call dealloc_sparse(mat_csr_cmplx)
		!								   messages,filesave,loadsym,savenum
		        call umfpack_unsym(mat_csc_cmplx%ni, mat_csc_cmplx, rhs(:,cmode), uvp, .false.,1,.false., .false.)

		        call dealloc_sparse(mat_csc_cmplx)
		        deallocate(rhs)
		elseif (P%solver .eq. 'pardiso') then
		!	USE PARDISO TO SOLVE THE SYSTEM
!			  if ((P%pardiso_iterative == 1).and.(itn > 1)) CGS_ready = ( max( delta(itn-1, 7, ccpt), delta(itn-1, 8, ccpt) ) < 0.1)
			  CGS_ready = .false. ! do not use the iterative solver
			  iter = ccpt*cmode ! So that reorderring is done every time when switching from Baro to ITM solvers
		      call pardiso_unsym(mat_csr_cmplx%ni, mat_csr_cmplx, rhs(:,cmode), uvp, P, iter, &
		      									 CGS_ready, dir_mats // 'temp/') ! directory where description file (dir_global) and factors will be saved

			  ! mat_csr_cmplx and rhs are DEALLOCATED inside pardiso_unsym(...)
	!	      deallocate(rhs)
		elseif (P%solver .eq. 'gmres') then
			  iter = 2 ! only matters that iter > 1 and pardiso has been initialized
			  allocate(rhs_cur(np), stat=istat)
			  rhs_cur = rhs(:,cmode)
			  call mkl_gmres(mat_csr_cmplx%ni, mat_csr_cmplx, rhs_cur, uvp, P, iter, cpt, dir_sols)
			  ! mat_csr_cmplx and rhs are DEALLOCATED inside mkl_gmres(...)
		else
		      print *, "You must choose a valid solver: pardiso or umfpack"
		      stop
		end if
	!********************************************************************
			if (P%messages >= 1) then
				call CPU_Time(T2_cpt)
				call system_clock ( wall_t2_cpt, clock_rate, clock_max )
				call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2_cpt-T1_cpt))) //', Wall: '&
							//trim(ADJUSTL(conv_secs( real(wall_t2_cpt-wall_t1_cpt)/real(clock_rate) )))//')')
			endif
		!==========================================================================================
		!    %================
		!    % update solution
		!    %================
			call save_itm_solution(nu, nv, nh, cmode, cpt, uvp, P%p_itm_avg, dir_sols//'itm/')
			deallocate(uvp)
		enddo

	deallocate(issponge_u, issponge_v, issponge_h) ! reload for another constituent

enddo

!=========================================================================================
!=================================DONE====================================================
!=========================================================================================
!    %==============================
!    % plot global tides and errors
!    %==============================
!
if (1==0) then
         call load_alloc_matrix(hp, dir_cols // 'hp.dat')
         call load_alloc_matrix(up, dir_cols // 'up.dat')
         call load_alloc_matrix(vp, dir_cols // 'vp.dat')

	do ccpt = 1, ncpts

	     cpt=cpts(2*ccpt-1:2*ccpt)
	     do cmode = 1, nmodes
		     call load_alloc_vector(u, dir_sols//'itm/' // cpt // '_u' // '_m'//tostring(cmode) // '.dat')
		     call load_alloc_vector(v, dir_sols//'itm/' // cpt // '_v' // '_m'//tostring(cmode) // '.dat')
		     call load_alloc_vector(pr, dir_sols//'itm/' // cpt // '_p' // '_m'//tostring(cmode) // '.dat')

	     	write(*,*) cpt,' mode', cmode, " max internal tide pressure [N/m2]: ", P%rhoo*maxval(abs(pr))
	     	write(*,*) cpt,' mode', cmode, " max modal velocity amp [m/s]: ", max(maxval(abs(u)),maxval(abs(v)))

			call save_sol4plot_cmplx(nph, nth, nh, hp, pr, dir_sols //'itm/' // 'sol_p_'//cpt//'_m'//tostring(cmode) // '.dat')
			call save_sol4plot_cmplx(nph, nth, nu, up, u, dir_sols //'itm/' // 'sol_u_'//cpt//'_m'//tostring(cmode) // '.dat')
			call save_sol4plot_cmplx(nph, nth, nv, vp, v, dir_sols //'itm/' // 'sol_v_'//cpt//'_m'//tostring(cmode) // '.dat')
		enddo
	enddo

	deallocate(hp, up, vp)

endif

end subroutine itm_solver_linear
!==========================================================================================

!==========================================================================================

subroutine itm_rhs(rhs, cpt, P, GD, nmodes, dir_cols, dir_grid, dir_mats, dir_sols)

    implicit none
!
!% calculates the barotropic tidal forcing for IT generation equation

     character(len=2), intent(in)   :: cpt
     type(params), intent(in)		:: P
     type(grid_dims), intent(in)	:: GD
     integer, intent(in)			:: nmodes
     integer 						:: nph, nth, nu, nv, nh


     integer  :: cooruv
     real(dp), allocatable		:: dc2n(:,:)
	 real(wp), allocatable		:: issponge_h(:, :)
     real(wp), allocatable		:: gradH_u(:), gradH_v(:), H_u(:), H_v(:), ta_u(:), ta_v(:)
     type (triplet) :: u2h, v2h

     complex(cwp), allocatable  :: tempu(:), tempv(:), tmp1(:), tmp2(:)

     integer ::      istatus, j, nf, ch, istat
!     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     complex(cwp), allocatable :: u0(:), v0(:)
     complex(cwp), allocatable :: rhs(:,:)

     real(wp), allocatable		:: metrics_hg(:, :)
	 integer, allocatable    	:: hp(:, :), mask_hg(:,:), mask_hg_sponge(:,:)
     complex(cwp), allocatable	:: rhs_hg(:,:)

     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols


	nph = GD%nph
	nth = GD%nta
	nu = GD%nu
	nv = GD%nv
	nh = GD%nh

!%====================================================================
!% differentiate between lat-lon/MERC, velocity/transport formulations
!%====================================================================
!	Only Transport-uvh formulation so far
if (P%coor == 1) then ! MERCATOR
!  if (P%uvform == 1) then
!    cooruv=11 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=12 ! TRANSPORT
!  endif
elseif (P%coor == 2) then ! LAT/LON
!  if (P%uvform == 1) then
!    cooruv=21 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=22 ! TRANSPORT
!  endif
endif

!%==============================================================
!% Initialize rhs
!%==============================================================
     allocate(rhs(nu+nv+nh, nmodes), stat = istatus)
     rhs = 0.
	! and auxillary things
	allocate(tempu(nu), tempv(nv),tmp1(nh), tmp2(nh), stat = istatus)
!%==================
!% Load grad(H) (calculated in write_baro_mats)
!%==================
	call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
	call load_alloc_vector(u0, dir_sols // cpt // '_u' // '_m0' // '.dat')
	call load_alloc_vector(v0, dir_sols // cpt // '_v' // '_m0' // '.dat')

	!	convert transport back to U
	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
			call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')
			tempu=(u0 * cosh(ta_u) / H_u)
			tempv=(v0 * cosh(ta_v) / H_v)
			deallocate(ta_u, ta_v)
		case (22) ! LAT-LON, TRANSPORT
			tempu=(u0 / H_u)
			tempv=(v0 / H_v)
	end select
	deallocate(u0, v0, H_u, H_v)

	call load_alloc_vector(gradH_u, dir_cols // 'gradH_u.dat')
	call load_alloc_vector(gradH_v, dir_cols // 'gradH_v.dat')
	! Load the necessary matrices: v2h, u2h
	call load_alloc_sparse(v2h, dir_mats // 'v2h.dat')
	call load_alloc_sparse(u2h, dir_mats // 'u2h.dat')

	! Calculate the dot product (u0, grad(H))
	call coo_vec_mul(u2h, tempu*cmplx(gradH_u, 0._dp, kind=cwp), P%lib, tmp1)
	call coo_vec_mul(v2h, tempv*cmplx(gradH_v, 0._dp, kind=cwp), P%lib, tmp2)

	deallocate(gradH_u, gradH_v)
	call dealloc_sparse(v2h)
	call dealloc_sparse(u2h)

	!***************************************************
	nf=ceiling(.5*(P%itm_avg*P%nph/360 - 1.)) ! will average over 2*nf+1 cells (extra nf in ph and th directions)
	! calculate metrics and masks to later apply running average (u- and v- grids)
	if ( (P%smooth_forcing==1).and.(nf >= 1) ) then
		if (P%messages >= 1) then
			write(*, '("(smoothed over a ", a, " grid pts per edge square) ")', advance='no') tostring(2*nf+1)
		endif
	    call load_alloc_matrix(hp, dir_cols // 'hp.dat')

	    call load_alloc_matrix(metrics_hg, dir_grid // 'temp/metrics_ug.dat')
		call load_alloc_matrix(mask_hg, dir_grid // 'temp/mask_hg.dat')
		allocate (rhs_hg(nph, nth), stat=istat)
	    if (P%sponge_forcing == 0) then
	    	allocate (mask_hg_sponge(nph, nth), stat=istat)
	    endif
	 endif
	!***************************************************

	! now calc rhs as: −(c2n)' ( U0 · grad(H) )

	call load_alloc_matrix(dc2n, dir_grid// 'itm/' // 'dc2n_h.dat')
	call load_alloc_matrix(issponge_h, dir_grid//'itm/'//'issponge_h'//cpt(2:2)//'.dat')

	do j = 1, nmodes

		rhs(nu+nv+1:nu+nv+nh, j) = -dc2n(:, j)*(tmp1 + tmp2)
		if (P%sponge_forcing == 0) then
			where (issponge_h(:,j)>0)	rhs(nu+nv+1:nu+nv+nh, j) = 0
		endif

		! Save original rhs before smoothing
!		call save_vector(rhs(nu+nv+1:nu+nv+nh, 1), dir_cols//'itm/' //'rhs_'// cpt // '_orig'// '.dat')

		if ( (P%smooth_forcing==1).and.(nf >= 1) ) then

			rhs_hg = 0
			do ch = 1, nh
			    rhs_hg(hp(ch,1),hp(ch,2)) = rhs(nu+nv+ch, j)
			enddo

			if (P%sponge_forcing == 0) then
				mask_hg_sponge = 0
				do ch = 1, nh
					mask_hg_sponge(hp(ch,1),hp(ch,2)) = log2int(issponge_h(ch,1)==0)
				enddo
			endif

		! Apply running average to smooth baro forcing (necessary due to nasty small-scale param ITD)
			if (P%sponge_forcing == 0) then
				call mv_average(rhs_hg, mask_hg*mask_hg_sponge, metrics_hg, nf, P%coor, P%lib)
			else
				call mv_average(rhs_hg, mask_hg, metrics_hg, nf, P%coor, P%lib)
			endif
		! index back onto an array
			do ch = 1, nh
			    rhs(nu+nv+ch, j) = rhs_hg(hp(ch,1),hp(ch,2))
			enddo
		endif
	enddo

	deallocate(dc2n,issponge_h)
!     call save_matrix(rhs(nu+nv+1:nu+nv+nh, :), dir_sols //'itm/'// 'rhs_'//cpt//'.dat')

end subroutine itm_rhs

!**********************************************************************************************************
!**********************************************************************************************************
subroutine itm_drag_m0(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

!% writes to file the drag of nmodes for all the baro constituents

implicit none

!     type(params), intent(in) :: P
!     integer            :: cooruv
     type(params), intent(in)	:: P
     type(grid_dims), intent(in):: GD
     character(len=*), intent(in)   :: cpts
     character(len=2)   :: cpt

     integer :: nph, nth, nu, nv, nh, nmodes, nf
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max
     real    ::      T1_cpt, T2_cpt ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1_cpt, wall_t2_cpt

     integer            :: istat, ccpt, cmode, ncpts, cu,cv, j

     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols

     type (triplet) :: h2v, h2u
     real(wp), allocatable		:: gradH_u(:), gradH_v(:)
     complex(cwp), allocatable	:: pr(:), prstore_u(:), prstore_v(:), prstore_h(:)
	 real(wp), allocatable		:: issponge_h(:, :), issponge_u(:, :), issponge_v(:, :)
     complex(cwp), allocatable	:: Du(:), Dv(:)

!     real(wp), allocatable		:: H_u(:), H_v(:), H_ug(:, :), H_vg(:, :),
!	 real(wp), allocatable		:: th_ug(:), ta_ug(:), th_vg(:), ta_vg(:)
     real(wp), allocatable		:: metrics_ug(:, :), metrics_vg(:, :)
	 integer, allocatable    	:: up(:, :), vp(:, :), mask_ug(:,:), mask_vg(:,:), mask_ug_sponge(:,:), mask_vg_sponge(:,:)
	complex(cwp), allocatable	:: u(:), v(:)
     complex(cwp), allocatable	:: ug(:,:), vg(:,:), Dug(:,:), Dvg(:,:)
     complex(cwp), allocatable	:: beta_u(:), beta_v(:), beta_ug(:,:), beta_vg(:,:)

! shortcuts
	nmodes = P%n_modes

	nph = GD%nph
	nth = GD%nta
	nu = GD%nu
	nv = GD%nv
	nh = GD%nh

	ncpts=len(cpts)/2
!%=======================
!% calculate the Q and D
!%=======================
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

write(*, '("  Storing IT drag: ")', advance='no')

	call load_alloc_vector(gradH_u, dir_cols // 'gradH_u.dat')
	call load_alloc_vector(gradH_v, dir_cols // 'gradH_v.dat')
	! Load the necessary matrices: h2u, h2v
	call load_alloc_sparse(h2u, dir_mats // 'h2u.dat')
	call load_alloc_sparse(h2v, dir_mats // 'h2v.dat')

	!***************************************************
	nf=ceiling(.5*(P%itm_avg*P%nph/360 - 1.)) ! will average over 2*nf+1 cells (extra nf in ph and th directions)
	! calculate metrics and masks to later apply running average (u- and v- grids)
	if ( (modulo(P%itm_method,10)>=1).and.(nf >= 1) ) then
		if (P%messages >= 1) then
			write(*, '("(smoothed over a ", a, " grid pts per edge square) ")', advance='no') tostring(2*nf+1)
		endif
	    call load_alloc_matrix(up, dir_cols // 'up.dat')
	    call load_alloc_matrix(vp, dir_cols // 'vp.dat')

	    call load_alloc_matrix(metrics_ug, dir_grid // 'temp/metrics_ug.dat')
		call load_alloc_matrix(metrics_vg, dir_grid // 'temp/metrics_vg.dat')
		call load_alloc_matrix(mask_ug, dir_grid // 'temp/mask_ug.dat')
		call load_alloc_matrix(mask_vg, dir_grid // 'temp/mask_vg.dat')
!
		allocate (Dug(nph, nth),Dvg(nph, nth+1), stat=istat)
		if ( (modulo(P%itm_method,10)==1)) then
			allocate (ug(nph, nth),vg(nph, nth+1), stat=istat)
		endif

	    if (P%sponge_forcing == 0) then
	    	allocate (mask_ug_sponge(nph, nth), mask_vg_sponge(nph, nth+1), stat=istat)
	    endif

	 endif
	!***************************************************

allocate(prstore_u(nu),prstore_v(nv),prstore_h(nh), stat=istat)
allocate(Du(nu), Dv(nv), stat = istat)
if ( (modulo(P%itm_method,10)==1).and.(nf >= 1) ) then
	allocate(beta_u(nu), beta_v(nv), stat = istat)
endif

do ccpt = 1, ncpts
     cpt=cpts(2*ccpt-1:2*ccpt)

		if (P%messages >= 1) then
			write(*, '(a, ": ")', advance='no') cpt
		endif

	if (P%sponge_forcing == 0) then
!		if (allocated(issponge_h)) deallocate(issponge_h)
		call load_alloc_matrix(issponge_h, dir_grid//'itm/'//'issponge_h'//cpt(2:2)//'.dat')
	endif

	 prstore_h = 0
     do cmode = 1, nmodes
	     call load_alloc_vector(pr, dir_sols//'itm/' // cpt // '_p' // '_m'//tostring(cmode) // '.dat')
		if (P%sponge_forcing == 0) then
			where (issponge_h(:,cmode)/=0)
				pr = 0
			end where
		endif

		prstore_h = prstore_h + pr
	 enddo! Note: multiply prstore by rho to get actual baroclinic pressure

	! Calculate pr_store on u and v grids
	call coo_vec_mul(h2u, prstore_h, P%lib, prstore_u)
	call coo_vec_mul(h2v, prstore_h, P%lib, prstore_v)
	! Now mask out

	! Calculate IT drag
	Du = -gradH_u*prstore_u
	Dv = -gradH_v*prstore_v

	if ( (modulo(P%itm_method,10)>0).and.(nf >= 1) ) then
		!***************************************************
		! Save original Du, Dv before smoothing
!		call save_vector(Du, dir_sols//'itm/' // cpt // '_Du' // '_m0_itd_drag_orig'// '.dat')
!		call save_vector(Dv, dir_sols//'itm/' // cpt // '_Dv' // '_m0_itd_drag_orig'// '.dat')

		!***************************************************
		!	PERFORMANCE: START
		!***************************************************
		if (P%messages >= 2) then
			write(*, '(a)', advance='no') 'moving average... '
		endif
		call CPU_Time(T1_cpt)
		call system_clock ( wall_t1_cpt, clock_rate, clock_max )

		Dug = 0; Dvg = 0;
		do cu = 1, nu
		    Dug(up(cu,1),up(cu,2)) = Du(cu)
		enddo
		do cv = 1, nv
		    Dvg(vp(cv,1),vp(cv,2)) = Dv(cv)
		enddo

		if (P%sponge_forcing == 0) then
!			if (allocated(issponge_u)) deallocate(issponge_u)
			call load_alloc_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//cpt(2:2)//'.dat')
			mask_ug_sponge = 0
			do cu = 1, nu
				mask_ug_sponge(up(cu,1),up(cu,2)) = log2int(issponge_u(cu,1)==0)
			enddo
!			if (allocated(issponge_v)) deallocate(issponge_v)
			call load_alloc_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//cpt(2:2)//'.dat')
			mask_vg_sponge = 0
			do cv = 1, nv
				mask_vg_sponge(vp(cv,1),vp(cv,2)) = log2int(issponge_v(cv,1)==0)
			enddo
		endif

			! calculate "beta" for IT drag
		if ( (modulo(P%itm_method,10)==1)) then
			!***************************************************
			call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
			call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')

			ug = 0; vg = 0;
			do cu = 1, nu
			    ug(up(cu,1),up(cu,2)) = u(cu)
			enddo
			do cv = 1, nv
			    vg(vp(cv,1),vp(cv,2)) = v(cv)
			enddo
	!print *, "hi3"
			if (P%sponge_forcing == 0) then
				call calc_sm_beta(Dug, ug, mask_ug*mask_ug_sponge, metrics_ug, nf, beta_ug, P%lib)
				call calc_sm_beta(Dvg, vg, mask_vg*mask_vg_sponge, metrics_vg, nf, beta_vg, P%lib)
			else
				call calc_sm_beta(Dug, ug, mask_ug, metrics_ug, nf, beta_ug, P%lib)
				call calc_sm_beta(Dvg, vg, mask_vg, metrics_vg, nf, beta_vg, P%lib)
			endif

			do cu = 1, nu
			    beta_u(cu) = beta_ug(up(cu,1),up(cu,2))
			enddo
			do cv = 1, nv
			    beta_v(cv) = beta_vg(vp(cv,1),vp(cv,2))
			enddo

!			if (P%sponge_forcing == 0) then
!				if (allocated(issponge_u)) deallocate(issponge_u)
!				call load_alloc_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//cpt(2:2)//'.dat')
!				if (allocated(issponge_v)) deallocate(issponge_v)
!				call load_alloc_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//cpt(2:2)//'.dat')
!
!				where (issponge_u(:,1)/=0) beta_u = 0 ! first mode gives the smallest sponge
!				where (issponge_v(:,1)/=0) beta_v = 0
!			endif
			!***************************************************
			call save_vector(beta_u, dir_sols//'itm/temp/' // cpt // '_beta_u' // '_m0_itd_drag'// '.dat')
			call save_vector(beta_v, dir_sols//'itm/temp/' // cpt // '_beta_v' // '_m0_itd_drag'// '.dat')

		elseif (modulo(P%itm_method,10)==2) then
		! Apply running average to smooth IT drag
			if (P%sponge_forcing == 0) then
!				print *, size(mask_ug),size(mask_vg), size(mask_ug_sponge),size(mask_vg_sponge)
				call mv_average(Dug, mask_ug*mask_ug_sponge, metrics_ug, nf, P%coor, P%lib)
				call mv_average(Dvg, mask_vg*mask_vg_sponge, metrics_vg, nf, P%coor, P%lib)
			else
				call mv_average(Dug, mask_ug, metrics_ug, nf, P%coor, P%lib)
				call mv_average(Dvg, mask_vg, metrics_vg, nf, P%coor, P%lib)
			endif

			do cu = 1, nu
			    Du(cu) = Dug(up(cu,1),up(cu,2))
			enddo
			do cv = 1, nv
			    Dv(cv) = Dvg(vp(cv,1),vp(cv,2))
			enddo

!			if (P%sponge_forcing == 0) then
!				if (allocated(issponge_u)) deallocate(issponge_u)
!				call load_alloc_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//cpt(2:2)//'.dat')
!				if (allocated(issponge_v)) deallocate(issponge_v)
!				call load_alloc_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//cpt(2:2)//'.dat')
!
!				where (issponge_u(:,1)/=0) Du = 0 ! first mode gives the smallest sponge
!				where (issponge_v(:,1)/=0) Dv = 0
!			endif

		endif

!print *, "hi4"
		!***************************************************
		!	PERFORMANCE: END
		!***************************************************
		  call CPU_Time(T2_cpt)
		     call system_clock ( wall_t2_cpt, clock_rate, clock_max )

		if (P%messages >= 2) then
			call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2_cpt-T1_cpt))) //', Wall: '&
						//trim(ADJUSTL(conv_secs( real(wall_t2_cpt-wall_t1_cpt)/real(clock_rate) )))//')')
		elseif (P%messages == 0) then
			write(*, '("")')
		endif

	endif

	!***************************************************
	call save_vector(Du, dir_sols//'itm/' // cpt // '_Du' // '_m0_itd_drag'// '.dat')
	call save_vector(Dv, dir_sols//'itm/' // cpt // '_Dv' // '_m0_itd_drag'// '.dat')

	!***************************************************
	if (P%messages == 1) then
     write(*, '("*,")', advance='no')
    endif

enddo

if ( (modulo(P%itm_method,10)>0).and.(nf >= 1) ) then
	deallocate(up, vp, Dug, Dvg)
	if ( (modulo(P%itm_method,10)==1)) then
		deallocate(u, v)
		deallocate(ug, vg)
	endif
endif
deallocate(prstore_u, prstore_v, prstore_h)
deallocate(Du, Dv)

  call CPU_Time(T2)
  call system_clock ( wall_t2, clock_rate, clock_max )
	if (P%messages >= 1 ) then
		write(*, '(a, a, a, a, a)') ' done (CPU: ', trim(ADJUSTL(conv_secs(T2-T1))), ', Wall: ',&
						trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) ))), ')'
	elseif (P%messages == 0 ) then
		write(*, '(" ")')
	endif

end subroutine itm_drag_m0
!**********************************************************************************************************
!**********************************************************************************************************

end module itm_solver_mod



