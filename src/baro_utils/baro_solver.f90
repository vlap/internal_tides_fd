module baro_solver_mod

!     use my_trigs
!     use generate_global_matrices

	 use itm_solver_mod
	 use generate_grid
	 use control
	 use baro_integrals
     use baro_projection
     use unsym_solvers
     use iter_solvers
     use my_sparse_aggregate
     use my_sparse
     use sal
     use itd
     use save_load
     use dispmodule
     use precisions, only: wp, cwp

     !use f90_kind      ! Kind-definition MODULE of the F-compiler
!     use sparse_utils  ! module for sparse matrix operations

!       External Subroutines
!      external   dzasum
!      real       dzasum

     contains

!==========================================================================================
subroutine baro_solver(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols, matlab_dir, save_gridid) ! dir_global

    implicit none
!%   Calculates a barotropic tide (iterative). Uses u/v/h + volume transport formulation
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points
type(disp_settings) ds

     integer ::     nph, nth, coor
     integer :: nu, nv, nh, np
     real(wp) :: latP, lonP
     real(wp) :: g, re, cdg, Q
     real(wp), allocatable :: beta(:), beta0(:)
!     complex(cwp) :: beta0_cmplx
!     complex(cwp), allocatable :: beta_cmplx(:)
     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols, matlab_dir, save_gridid!, dir_global
     type (tide_params) :: pars
     type(params)	:: P
     type(grid_dims):: GD

	 logical :: 	 CGS_ready = .false.!, last_itn_flag = .false.
!	 logical :: 	 file_exist = .false.
	 integer ::      itn, itn_reorder, cvgd, itm_active
     integer ::      j, istat, l,m,n
     real    ::      T1_itn, T2_itn, T1_cpt, T2_cpt ! for measuring CPU (NOT REAL TIME!)
     integer :: 	 wall_t1_itn, wall_t2_itn, wall_t1_cpt, wall_t2_cpt, clock_rate, clock_max
     integer, allocatable :: hp(:, :), up(:, :), vp(:, :)

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt
     integer :: iu(GD%nu), iv(GD%nv), ih(GD%nh)

     complex(cwp), allocatable :: tmp_u(:), tmp_v(:)
     complex(cwp), allocatable :: hsal_res(:)
     type (triplet) :: GHU, GHV
	 type (triplet_cmplx) :: DUU, DVU, DUV, DVV
     type (triplet) :: u2v, v2u, u2v_r, v2u_r, u2v_i, v2u_i

	 type (csr) :: tmp_csr
	 type (csr_cmplx) :: tmp_csr_cmplx

     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr) :: mat_csr
     type (csr_cmplx) :: mat_csr_cmplx
     type (csc_cmplx) :: mat_csc_cmplx

     type (csr_cmplx) :: speye
     integer, allocatable :: bcdiag(:)
     complex(cwp), allocatable, dimension(:) :: rhs

     complex(cwp), allocatable :: uvh(:), u(:), v(:), h(:), hsal(:)

     complex(cwp), allocatable :: Du(:), Dv(:), Qu(:), Qv(:), tmp(:), betau(:), betav(:)
	 real(wp), allocatable :: KU(:), KV(:)
     real(wp), allocatable     :: H_u(:), H_v(:)!, H_h(:)
	 real(wp), allocatable		:: gradH_u(:), gradH_v(:)
     logical, allocatable, dimension(:,:) :: u_crit_lat,v_crit_lat

     real(wp), allocatable     :: delta(:, :, :) !delta(itn,type,ccpt)

     type(domain_integrals) :: di

!==========================================================================================
! DISP MODULE SETTINGS
!==========================================================================================
!call tostring_set(rfmt='F12.1')
!==========================================================================================
! SHORTCUTS
!==========================================================================================
latP = P%latP
lonP = P%lonP
g = P%g
re = P%re
cdg = P%cd
Q = P%Qbar

nph = GD%nph
nth = GD%nta
nu = GD%nu
nv = GD%nv
nh = GD%nh
np = GD%np
!==========================================================================================

write(*, '("====================================")')
write(*, '("Welcome to baro_solver (nonlinear)")')
write(*, '("====================================")')

itn=1         ! iteration number
cvgd=0        ! switched to 1 when convergence occurs
if ((P%fr_scheme <= 2).and.(P%sal_scheme <= 1).and.(P%itd_scheme >= 2)) then ! Baro problem is linear. Turn on IT modelling at itn = 1
	itm_active=1 ! % switched to 1 when local modelling starts
else
	itm_active=0 ! % 0 otherwise
endif

ncpts=len(cpts)/2

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (nu + j), j=1,nv ) /)
     ih = (/ ( (nu + nv + j), j=1,nh ) /)

if (P%sal_scheme >= 2) then
	allocate(beta0(ncpts), stat = istat)
	beta0 = P%beta0
endif

! precalculate u_crit_lat,v_crit_lat that indicate regions of radiating/coastally trapped ITs
! (u_crit_lat,v_crit_lat) are .TRUE. within the (-lat_critical, lat_critical) where ITs radiate away
if ( (P%itd_scheme >= 2).and.(P%itm_method/10==1) ) then
	call calc_crit_lat_prefactor(u_crit_lat,v_crit_lat, cpts, P, nu, nv, dir_cols, dir_grid)
endif

!

do while ( (cvgd == 0) .and. (itn < P%itn_max)) !.and. (maxval(delta(itn,3:6,:)) < 1000) )

	allocate(delta(itn, 8, len(cpts)/2), stat = istat)
	delta=0.
	if (itn > 1) then
	    call load_matrix_3d(delta(1:itn-1, :, :), dir_sols // 'delta.dat')
	endif

	  call disp('==============')
	  call disp(' Iteration ', itn)
	  call disp('==============')

	write(*, '(" ------------ BAROTROPIC TIDE ------------ ")')

	!%==========================
	!% solve for each component
	!%==========================
	call CPU_Time(T1_itn)
	call system_clock ( wall_t1_itn, clock_rate, clock_max )

	do ccpt = 1, ncpts

	    cpt=cpts(2*ccpt-1:2*ccpt)

	!%==================================================
	!% Load the basic sparse matrix generated by baro_uvhT_mat
	!% (includes terms: rotation, c.o.mass, itd and RHS (forcing))
	!%==================================================
	!     call load_alloc_cmat(mat, dir_mats // 'temp/mat_init.dat')
	     ! Convert mat from COO to CSR:
	!     call coo2csr(mat, mat_csr, P%lib)
	!     call dealloc_sparse(mat)
		if ( (P%itd_scheme == 0).or.((P%itd_scheme == 2) .and. (itm_active == 1)) ) then
			call load_alloc_sparse(mat_csr, dir_mats // 'temp/' // 'mat_csr_init.dat')
		else
			call load_alloc_sparse(mat_csr, dir_mats // 'temp/' // cpt // '_mat_csr_init.dat')
		endif
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

	!     mat = mat-i*omega0*speye(np)
	      pars = get_pars(cpt)
	      speye = init_speye(np, cmplx(0., -pars%omega0, kind=cwp) )
	      call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), speye, mkl)
	      call dealloc_sparse(speye)

	!  %=====================
	!  % Print matrix dimensions
	!  %=====================
	    if ((itn == 1).and.(ccpt == 1)) then
	    	call tostring_set(rfmt='ES0.2')
	      	call disp('Have [' // tostring(real(mat_csr_cmplx%ni)) //' x '// tostring(real(mat_csr_cmplx%nj))//'] sparse matrix ('// &
	                 'nnz = ' // tostring(real(mat_csr_cmplx%nz))//'). Solve with '//P%solver// '.' )
			call tostring_set_factory()
	    endif

	write(*, '(i2, ". for ", a, ": Assembling the matrix and rhs, ")', advance = 'no') ccpt, cpt
	!%=======================
	!% calculate rhs forcing
	!%=======================
		if (itn == 1) then
			call baro_rhs(rhs, cpt, P, GD, dir_cols, dir_grid, dir_mats)
			call save_vector(rhs, dir_cols //'temp/'// cpt // '_rhs_0.dat')
		else
			call load_alloc_vector(rhs, dir_cols //'temp/'// cpt // '_rhs_0.dat')
		endif
	!==========================================================================================
	!    %==============================================================
	!    % add friction  (already included if fr_scheme = 0 or 1 or 2)
	!    %==============================================================
	     if (P%fr_scheme >= 3) then

			call load_alloc_vector(KU, dir_cols // 'KU.dat')
			call load_alloc_vector(KV, dir_cols // 'KV.dat')
			allocate(Qu(nu), Qv(nv), stat = istat)

			if ( (itn == 1).and.(P%load_sol == 0)) then ! Use linear friction option 2 (average value for transport Q)
				Q = P%Qbar
				call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
				call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
				Qu = cdg*Q/H_u**2 / KU
				Qv = cdg*Q/H_v**2 / KV
				deallocate(H_u, H_v)
			else ! nonlinear bottom drag
			     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
			     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
			     call load_alloc_vector(Du, dir_sols // cpt // '_Du' // '_m0_drag' // '.dat')
			     call load_alloc_vector(Dv, dir_sols // cpt // '_Dv' // '_m0_drag' // '.dat')

			   where (u==0)
				   Qu=Du
			   elsewhere
				   Qu=Du/u
			    end where
			   where (v==0)
				   Qv=Dv
			   elsewhere
				   Qv=Dv/v
			    end where
	!		   Du=Du-P%p_fric*Qu*u
	!		   Dv=Dv-P%p_fric*Qv*v
	!		   Qu=P%p_fric*Qu
	!		   Qv=P%p_fric*Qv

			     rhs(iu) = rhs(iu) - KU*(Du - Qu*u)
			     rhs(iv) = rhs(iv) - KV*(Dv - Qv*v)

			   deallocate(u,v, Du, Dv)
			endif

	!       mat = mat + sparse(iu,iu,Qu./cosh(ta_u),np,np) + sparse(iv,iv,Qv./cosh(ta_v),np,np);
	!	    call init_sparse_vals(speye,(/ ( (j), j=1,nu+nv+1), ((nu+nv+1), j=nu+nv+1, np+1) /),(/ ( (j), j=1,nu+nv) /), &
	!		           (/ KU*Qu, KV*Qv/), np, np, nu+nv) ! Initialize diagonal CSR matrix
	!		      call dealloc_sparse(speye)
	!		call disp(speye%vals, orient='row')
			allocate(tmp(nh), stat=istat)
			tmp = 0.
			call csr_dia_add(mat_csr_cmplx, (/  KU*Qu, KV*Qv, tmp/))
			deallocate(tmp)

			deallocate(KU, KV, Qu, Qv)
		endif
	!==========================================================================================

	!==========================================================================================
	!    %==============================================================
	!    % add IT drag  (already included if itd_scheme = 0 or 1 )
	!    %==============================================================
	     if ( (P%itd_scheme >= 2).and.(itm_active==1).and.(itn > 1) ) then

			call load_alloc_vector(KU, dir_cols // 'KU.dat')
			call load_alloc_vector(KV, dir_cols // 'KV.dat')
			allocate(Qu(nu), Qv(nv), stat = istat)

			call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
			call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
			call load_alloc_vector(Du, dir_sols//'itm/' // cpt // '_Du' // '_m0_itd_drag'// '.dat')
			call load_alloc_vector(Dv, dir_sols//'itm/' // cpt // '_Dv' // '_m0_itd_drag'// '.dat')

		! choose the scheme with which dynamical IT drag is incorporated into the baro solve
			if (P%itm_method/10==1) then

!				inquire( file=dir_sols//'itm/temp/' // cpt // '_beta_u' // '_m0_itd_drag'// '.dat', exist=file_exist )
				if (P%itm_method-10==1) then
					call load_alloc_vector(Qu, dir_sols//'itm/temp/' // cpt // '_beta_u' // '_m0_itd_drag'// '.dat')
					call load_alloc_vector(Qv, dir_sols//'itm/temp/' // cpt // '_beta_v' // '_m0_itd_drag'// '.dat')
				else
					!*****************
					! IMPORTANT NOTE:
					!*****************
					! Taking Qu = real(Du) works best within the (-lat_critical, lat_critical) where ITs radiate away
					! using Qu=Du/u there creates small-scale oscillations around lines where abs(u) is small and dynamical IT drag has relatively strong influence
					! However, the picture is reversed for the coastally trapped waves (CTW) closer to the poles. For them Qu = real(Du) will never converge!
					! hence, Qu=Du/u is taken and it seems to work well.
					!********************************************************************
				   ! u
				   where (u == 0)
					   Qu = Du
				   elsewhere (u_crit_lat(ccpt,:))	! for radiating internal waves
					   Qu = real(Du/u)
				   elsewhere						! for coastally trapped waves
					   Qu = Du/u
				   end where
				   ! v
				   where (v == 0)
					   Qv = Dv
				   elsewhere (v_crit_lat(ccpt,:))	! for radiating internal waves
					   Qv = real(Dv/v)
				   elsewhere						! for coastally trapped waves
					   Qv = Dv/v
				   end where
				endif
!				print *, "Qu cut-off (% of cells):", real(count(abs(Qu)>P%itm_qmax))*100/nu, &
!						 "Qv cut-off (% of cells):", real(count(abs(Qv)>P%itm_qmax))*100/nv

!			   where (real(Qu)<0) Qu = 0 ! does not work. blow up in negative conversion regions

			   where ((abs(Qu)>P%itm_qmax).and.(u_crit_lat(ccpt,:))) Qu=P%itm_qmax * Qu/abs(Qu)
			   where ((abs(real(Qu))>P%itm_qmax/50).and.(.not.u_crit_lat(ccpt,:))) &
			   					  Qu=P%itm_qmax/50 * real(Qu)/abs(real(Qu)) + imag(Qu)

!			   where (real(Qv)<0) Qv = 0 ! does not work. blow up in negative conversion regions

			   where ((abs(Qv)>P%itm_qmax).and.(v_crit_lat(ccpt,:))) Qv=P%itm_qmax * Qv/abs(Qv)
			   where ((abs(real(Qv))>P%itm_qmax/50).and.(.not.v_crit_lat(ccpt,:))) &
			   					  Qv=P%itm_qmax/50 * real(Qv)/abs(real(Qv)) + imag(Qv)

!				print *, 'rhs original'
!				print *, 'u', maxval(abs(rhs(iu))), maxval(real(rhs(iu))), minval(real(rhs(iu))), sum(rhs(iu))
!				print *, 'v', maxval(abs(rhs(iv))), maxval(real(rhs(iv))), minval(real(rhs(iv))), sum(rhs(iv))
				 rhs(iu) = rhs(iu) - KU*(Du - Qu*u)
			     rhs(iv) = rhs(iv) - KV*(Dv - Qv*v)
!				print *, 'rhs modified'
!				print *, 'u', maxval(abs(rhs(iu))), maxval(real(rhs(iu))), minval(real(rhs(iu))), sum(rhs(iu))
!				print *, 'v', maxval(abs(rhs(iv))), maxval(real(rhs(iv))), minval(real(rhs(iv))), sum(rhs(iv))

				allocate(tmp(nh), stat=istat)
				tmp = 0.
				call csr_dia_add(mat_csr_cmplx, (/  KU*Qu, KV*Qv, tmp/))
				deallocate(tmp)

			elseif (P%itm_method/10==2) then

				call load_alloc_vector(gradH_u, dir_cols // 'gradH_u.dat')
				call load_alloc_vector(gradH_v, dir_cols // 'gradH_v.dat')
				! Load the necessary matrices: v2u, u2v
				call load_alloc_sparse(v2u, dir_mats // 'v2u.dat')
				call load_alloc_sparse(u2v, dir_mats // 'u2v.dat')

				allocate(tmp_u(nu), tmp_v(nv), stat=istat)
				call coo_vec_mul(u2v, gradH_u*u, P%lib, tmp_v)
				call coo_vec_mul(v2u, gradH_v*v, P%lib, tmp_u)

				tmp_u = (u*gradH_u + tmp_u)*gradH_u
				tmp_v = (v*gradH_v + tmp_v)*gradH_v

			   where (tmp_u==0)
				   Qu=Du
!				   Qu=real(Du)
			   elsewhere
				   Qu=Du/tmp_u
!				   Qu=real(Du/tmp_u)
			   end where
			   ! adjust P%itm_qmax by 10^4 due to dividing by grad(H)^2
			   where (abs(Qu)>10**4 * P%itm_qmax) Qu=10**4 * P%itm_qmax * Qu/abs(Qu)

			   where (tmp_v==0)
				   Qv=Dv
!				   Qv=real(Dv)
			   elsewhere
				   Qv=Dv/tmp_v
!				   Qv=real(Dv/( tmp_v )) !
			   end where
!				print *, "Qu cut-off (% of cells):", real(count(abs(Qu)>10**4 * P%itm_qmax))*100/nu, &
!						 "Qv cut-off (% of cells):", real(count(abs(Qv)>10**4 * P%itm_qmax))*100/nv

			   where (abs(Qv)>10**4 * P%itm_qmax) Qv=10**4 * P%itm_qmax * Qv/abs(Qv)
			!*********************************
			!	RHS
			!*********************************
				 rhs(iu) = rhs(iu) - KU*(Du - Qu*tmp_u)
			     rhs(iv) = rhs(iv) - KV*(Dv - Qv*tmp_v)

			     deallocate(tmp_u, tmp_v)
			!*********************************
			!	Matrix
			!*********************************
			!	DUU
				call init_sparse_vals(DUU, iu, iu, KU*(Qu*gradH_u**2), np, np, nu)! an (np x np) matrix.
				! Convert the matrix from COO to CSR:
			    call coo2csr(DUU, tmp_csr_cmplx, P%lib)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr_cmplx)
			!	DVV
				call init_sparse_vals(DVV, iv, iv, KV*(Qv*gradH_v**2), np, np, nv)! an (np x np) matrix.
				! Convert the matrix from COO to CSR:
			    call coo2csr(DVV, tmp_csr_cmplx, P%lib)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr_cmplx)

			!	DVU
				call right_mult_diag(v2u, gradH_v, nv)
				call init_sparse_vals(v2u_r, v2u%indi, v2u%indj, v2u%vals, v2u%ni, v2u%nj, v2u%nz)
				call init_sparse_vals(v2u_i, v2u%indi, v2u%indj, v2u%vals, v2u%ni, v2u%nj, v2u%nz)

			    call left_mult_diag(v2u_r, KU*(real(Qu)*gradH_u), nu)
			    call left_mult_diag(v2u_i, KU*(imag(Qu)*gradH_u), nu)

				call init_sparse_vals(DVU, v2u%indi, v2u%indj+nu, cmplx(v2u_r%vals, v2u_i%vals, kind=cwp), np, np, v2u%nz)! an (np x np) matrix.
				! Convert the matrix from COO to CSR:
			    call coo2csr(DVU, tmp_csr_cmplx, P%lib)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr_cmplx)

			!	DUV
				call right_mult_diag(u2v, gradH_u, nu)
				call init_sparse_vals(u2v_r, u2v%indi, u2v%indj, u2v%vals, u2v%ni, u2v%nj, u2v%nz)
				call init_sparse_vals(u2v_i, u2v%indi, u2v%indj, u2v%vals, u2v%ni, u2v%nj, u2v%nz)

			    call left_mult_diag(u2v_r, KV*(real(Qv)*gradH_v), nv)
			    call left_mult_diag(u2v_i, KV*(imag(Qv)*gradH_v), nv)

				call init_sparse_vals(DUV, u2v%indi+nu, u2v%indj, cmplx(u2v_r%vals, u2v_i%vals, kind=cwp), np, np, u2v%nz)! an (np x np) matrix.
				! Convert the matrix from COO to CSR:
			    call coo2csr(DUV, tmp_csr_cmplx, P%lib)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr_cmplx)

			    ! Test:
!	call save_sparse(DUU, dir_mats // 'tUU.dat')
!	call save_sparse(DUV, dir_mats // 'tUV.dat')
!	call save_sparse(DVU, dir_mats // 'tVU.dat')
!	call save_sparse(DVV, dir_mats // 'tVV.dat')
				if (1==0) then ! do the test
				    allocate(tmp_u(np), tmp_v(np), stat=istat)

				    allocate(tmp(nh), stat=istat)
					tmp = 0.
				    call coo_vec_mul(DUU, (/u,v,tmp/), P%lib, tmp_u)
				    call coo_vec_mul(DVV, (/u,v,tmp/), P%lib, tmp_v)

	!			    print *, maxval(abs(Du)), maxval(abs(Dv))
	!			    print *, sum(abs(Du))/nu, sum(abs(Dv))/nv
				    Du = Du-tmp_u(1:nu)
				    Dv = Dv-tmp_v(nu+1:nu+nv)
	!			    print *, maxval(abs(Du)), maxval(abs(Dv))
	!			    print *, sum(abs(Du))/nu, sum(abs(Dv))/nv

				    call coo_vec_mul(DVU, (/u,v,tmp/), P%lib, tmp_u)
				    call coo_vec_mul(DUV, (/u,v,tmp/), P%lib, tmp_v)
				    Du = Du-tmp_u(1:nu)
				    Dv = Dv-tmp_v(nu+1:nu+nv)

	!			    print *, maxval(abs(Du)), maxval(abs(Dv))
	!			    print *, sum(abs(Du))/nu, sum(abs(Dv))/nv

				    deallocate(tmp, tmp_u, tmp_v)
			    endif

				call dealloc_sparse(u2v_r)
				call dealloc_sparse(v2u_r)
				call dealloc_sparse(u2v_i)
				call dealloc_sparse(v2u_i)
	!******************************************************************
				call dealloc_sparse(DUU)
				call dealloc_sparse(DUV)
				call dealloc_sparse(DVU)
				call dealloc_sparse(DVV)

				call dealloc_sparse(u2v)
				call dealloc_sparse(v2u)

				deallocate(gradH_u, gradH_v)

			endif

		    deallocate(u,v, Du, Dv)
			deallocate(KU, KV, Qu, Qv)
		endif
	!==========================================================================================

	!================================================================================================
	!    % add self-attraction and loading (already included if sal_scheme = 0 or 1)
	!================================================================================================
		if (P%sal_scheme >= 2) then
			!	Define: beta and hsal_res
			allocate(hsal_res(nh), beta(nh), stat = istat)

			if ((itn == 1).and.(P%load_sol==0)) then
				hsal_res = 0
				beta = P%beta0
			else
				call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')
				call load_vector(hsal_res, dir_sols //'temp/'// cpt // '_hsal' // '.dat') ! load hsal vector into hsal_res
				if (P%sal_scheme == 3) then
					call load_vector(beta, dir_sols //'temp/'// cpt // '_betasal' // '.dat')
				elseif (P%sal_scheme == 2) then
					beta = beta0(ccpt)
				endif

				hsal_res = hsal_res - beta * h

				deallocate(h)
		    endif
		!******************************************************************
		!		Add the corresponding terms into the linear matrix and rhs
		!		UVH-formulation only
		!******************************************************************
			!*********************************
			!	RHS
			!*********************************
			call load_alloc_sparse(GHU, dir_mats // 'GHU.dat')
			call load_alloc_sparse(GHV, dir_mats // 'GHV.dat')

			if (.not.((itn == 1).and.(P%load_sol==0))) then
	      		allocate(tmp_u(nu), tmp_v(nv), stat = istat)

				call coo_vec_mul(GHU, hsal_res, P%lib, tmp_u)
				rhs(iu) = rhs(iu) + tmp_u
				call coo_vec_mul(GHV, hsal_res, P%lib, tmp_v)
				rhs(iv) = rhs(iv) + tmp_v

				deallocate(tmp_u, tmp_v)
			end if
			!*********************************
			!	Matrix
			!*********************************
			!	pressure gradients: d/dph terms in the matrix
			    call right_mult_diag(GHU, - beta, nh)
			    ! Transform GHU into a (np x np) matrix.
			    GHU%ni = np
			    GHU%nj = np
			    GHU%indj = GHU%indj + nu + nv ! (move elements to the "h" columns of the matrix)

				! Convert the matrix from COO to CSR:
			    call coo2csr(GHU, tmp_csr, P%lib)
			    call init_sparse_vals(tmp_csr_cmplx,tmp_csr%indi,tmp_csr%indj,tmp_csr%vals, &
				                           tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr)
			    call dealloc_sparse(tmp_csr_cmplx)

			!	pressure gradients: d/dta terms in the matrix
			    call right_mult_diag(GHV, - beta, nh)
			    ! Transform GHV into a (np x np) matrix.
			    GHV%ni = np
			    GHV%nj = np
			    GHV%indi = GHV%indi + nu ! (move elements to the "v" rows of the matrix)
			    GHV%indj = GHV%indj + nu + nv ! (move elements to the "h" columns of the matrix)

				! Convert the matrix from COO to CSR:
			    call coo2csr(GHV, tmp_csr, P%lib)
			    call init_sparse_vals(tmp_csr_cmplx,tmp_csr%indi,tmp_csr%indj,tmp_csr%vals, &
				                           tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr)
			    call dealloc_sparse(tmp_csr_cmplx)
	!******************************************************************

		    call dealloc_sparse(GHU)
		    call dealloc_sparse(GHV)

			deallocate(hsal_res, beta)

		end if

	!==========================================================================================

	!==========================================================================================
	!    %================================================
	!    % add something to do with Coriolis matrices...
	!    %================================================

	!    if (f) > 1
	!      if itn > 1
	!        load([flags.dir.file,'/',domain,'/',cpt],'u','v');
	!        load(matfile,'v2uf','v2ufsp','u2vf','u2vfsp');
	!        rhs(iu)=rhs(iu)+(v2uf-v2ufsp)*v;
	!        rhs(iv)=rhs(iv)-(u2vf-u2vfsp)*u;
	!        clear u v v2uf v2ufsp
	!      end
	!    end

	!  %=============
	!  % impose bcs:
	!  %=============
	write(*, '("implementing b.c. ")')!, advance = 'no')

	call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
	     rhs = rhs * bcdiag
	!call save_vector(bcdiag, dir_cols // 'bcdiag' //'.dat')
	     deallocate(bcdiag)

	!write(*, '("DONE")')
	!==========================================================================================
	!     Save the final version of the matrix
	if (P%messages>=3) then
	      call csr2coo(mat_csr_cmplx, mat_coo_cmplx, P%lib)
	      call save_sparse(mat_coo_cmplx, dir_mats // 'mat_' // cpt //'.dat')
	      call dealloc_sparse(mat_coo_cmplx)
	endif

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
	        call umfpack_unsym(mat_csc_cmplx%ni, mat_csc_cmplx, rhs, uvh, .false.,1,.false., .false.)

	        call dealloc_sparse(mat_csc_cmplx)
	        deallocate(rhs)
	elseif (P%solver .eq. 'pardiso') then
	!	USE PARDISO TO SOLVE THE SYSTEM
		  if ((P%pardiso_iterative == 1).and.(itn > 1)) CGS_ready = ( max( delta(itn-1, 7, ccpt), delta(itn-1, 8, ccpt) ) < 0.1)

		  if ((P%itd_scheme >= 2) .and. (itm_active == 1)) then
				itn_reorder = ccpt ! So that reorderring is done every time when switching from Baro to ITM solvers
!				last_itn_flag = (ccpt==ncpts)
		  else
				itn_reorder = itn*ccpt ! Still solving for the barotropic tide only, reordering just for the first call of pardiso
!				last_itn_flag = .false.
		  endif
	     call pardiso_unsym(mat_csr_cmplx%ni, mat_csr_cmplx, rhs, uvh, P, itn_reorder, &
	      													 CGS_ready, dir_mats // 'temp/') ! directory where description file (dir_global) and factors will be saved

	      ! mat_csr_cmplx and rhs are DEALLOCATED inside pardiso_unsym(...)
!	      deallocate(rhs)
	elseif (P%solver .eq. 'gmres') then
		  call mkl_gmres(mat_csr_cmplx%ni, mat_csr_cmplx, rhs, uvh, P, itn, cpt, dir_sols)
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
		call save_solution(nu, nv, nh, itn, P%load_sol, cpt, uvh, P%p_avg, dir_sols)
		deallocate(uvh)
	!==========================================================================================
	!****************** CLEAN UP Memory and Hard Drive (Pardiso) ********
!	if  (P%solver .eq. 'pardiso') then
!		call pardiso_clean(mat_csr_cmplx%ni, P,  dir_mats // 'temp/')! directory where description file (dir_global) and factors will be saved
!	end if

	enddo
!===============================================================================================
	!  %==============================================
	!  % examine height differences, i.e. convergence
	!  %==============================================
	call calc_convergence(nu, nv, nh, itn, P, cpts, delta, dir_cols, dir_sols)

	  if (itn > 1) then
	    if ( (P%cvg_scheme == 1) .and. (maxval( delta(itn-1:itn,4,:) ) < P%cvg_dhbar) .and. &
	    								(maxval( delta(itn-1:itn,3,:) ) < P%cvg_dhmax) ) then
	        cvgd = 1
	    elseif ( (P%cvg_scheme == 2) .and. (maxval( delta(itn-1:itn,6,:) ) < P%cvg_dhbar) .and. &
	    									(maxval(delta(itn-1:itn,5,:) ) <P%cvg_dhmax) ) then
	        cvgd = 1
	    endif
	  endif
!==========================================================================================

!if (cvgd == 0) then
	!    %================================
	!    % calculating updated h_SAL terms
	!    %================================
	if (P%sal_scheme >= 2) then
		call calc_sal(ncpts, cpts, P, GD, dir_grid, dir_cols, dir_sols, beta0)
	endif

	!    %====================================
	!    % calculate updated nonlinear terms:
	!    %====================================
	if (P%fr_scheme >= 3) then
		call calc_dhat(cpts, nu, nv, P, dir_cols, dir_mats, dir_sols)
	endif
	!==========================================================================================

	!    %==========================
	!    % consider IT modelling ! The convergence tolerance for ITM modelling is 2^3 times the convergence tolerance
	!    %==========================
	if ((P%itd_scheme >= 2) .and. (itm_active == 0)) then
      	if ( (P%cvg_scheme == 1).and.(maxval( delta(itn,4,:) ) < 8*P%cvg_dhbar).and.(maxval( delta(itn,3,:) ) < 8*P%cvg_dhmax)) then
	        itm_active=1
	    elseif ((P%cvg_scheme == 2).and.(maxval( delta(itn,6,:) ) < 8*P%cvg_dhbar).and.(maxval(delta(itn,5,:)) < 8*P%cvg_dhmax)) then
	        itm_active=1
	    endif
	    if ( (itm_active==1).and.(P%baro_sol == 1) ) then
	    	! make a copy of the baro-only solution
			do ccpt = 1, ncpts
				cpt=cpts(2*ccpt-1:2*ccpt)
	    		call copy_baro_solution(cpt, dir_sols)
	    	enddo
	    endif

	    if ( (itm_active==1).and.(P%messages >= 1) ) then
!		calculate domain integrals ONCE when baro tide converged but before going into the IT modelling
			do ccpt = 1, ncpts
				cpt=cpts(2*ccpt-1:2*ccpt)
				call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, 0, dir_grid, dir_cols, dir_mats, dir_sols, di)
				call show_domain_integrals(cpt, di)
			enddo
!			P%messages = 2 ! temprorary to monitor IT modelling impact
		endif
    endif
! 	 %==============================
!    % IT modelling:
!    %==============================
	if ((P%itd_scheme >= 2) .and. (itm_active == 1)) then

		call forced_itm_solver(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

	endif

! 	 %=================================
!    % extra domain integrals each itn:
!    %=================================
	if ( (P%messages >= 2 ).or.(P%di_log==1) ) then
		do ccpt = 1, ncpts
			cpt=cpts(2*ccpt-1:2*ccpt)
		!	Print the integrals on every iteration
			call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, itm_active, dir_grid, dir_cols, dir_mats, dir_sols, di)
			call show_domain_integrals(cpt, di, P%di_log)
			if (P%di_log==1) then
				call log_di(cpt, itn, di, '/../../di_log.txt', dir_sols)
			endif
		enddo

	call disp('---------------------------------------------------------------------------')
	endif

!==========================================================================================

	call CPU_Time(T2_itn)
	call system_clock ( wall_t2_itn, clock_rate, clock_max )
	call disp ('===>>> Total time spent on Iteration ' // tostring(itn) //': ' &
	            // trim(conv_secs(T2_itn-T1_itn)) // ' CPU, ' &
	            //trim(conv_secs( real(wall_t2_itn-wall_t1_itn)/real(clock_rate) ))//'  Wall')
	!==========================================================================================

	!      %==============================
	!      % plot global tides and errors
	!      %==============================

		if  ( P%graphics == 2 ) then
	!		Use MATLAB functions to export the data in M-files into binaries
			call system('xterm -e matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
			'"addpath(genpath(''' // matlab_dir // ''')); baro_v1_plot_conv(''' // &
			trim(save_gridid)//''', '// tostring(itn) // ', ''' // cpts // ''', '//tostring(P%load_sol)//'); waitforbuttonpress; exit;" ')
	!		print *, ' matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
	!		'"addpath(genpath(''' // matlab_dir // ''')); save_grid_binary(''2012_06_13__16_54'')" '

	!		if the matlab script failed

	!		if  (.not. dir_e ) then
	!			write(*,'(a)') 'Directory /'//trim(P%gridid)//' for the previously generated grid files doesn''t exist.'
	!			stop
	!		end if
		end if
	!==========================================================================================

itn = itn + 1
!==========================================================================================
!endif


enddo

!=========================================================================================
!=================================DONE====================================================
!=========================================================================================
write(*, '("========================================================================")')
write(*, '("========================================================================")')


!    %==============================
!    % plot global tides and errors
!    %==============================
!
!do ccpt = 1, ncpts
!
!     cpt=cpts(2*ccpt-1:2*ccpt)
!
!     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
!     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
!     call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')
!
!     write(*,*) cpt, " max tide height: ", maxval(abs(h))
!
!         call load_alloc_matrix(hp, dir_cols // 'hp.dat')
!         call load_alloc_matrix(up, dir_cols // 'up.dat')
!         call load_alloc_matrix(vp, dir_cols // 'vp.dat')
!
!		call save_sol4plot_cmplx(nph, nth, nh, hp, h, dir_sols // 'sol_h_'//cpt//'.dat')
!		call save_sol4plot_cmplx(nph, nth, nu, up, u, dir_sols // 'sol_u_'//cpt//'.dat')
!		call save_sol4plot_cmplx(nph, nth, nv, vp, v, dir_sols // 'sol_v_'//cpt//'.dat')
!
!		deallocate(hp, up, vp)
!
!enddo

end subroutine baro_solver
!==========================================================================================

!==========================================================================================
subroutine baro_solver_linear(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

    implicit none
!%   Calculates a barotropic tide (linear bottom friction). Uses u/v/h + volume transport formulation
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points

     integer :: nph, nth
     integer :: nu, nv, nh, np
     real(wp):: latP, lonP
     real(wp):: g, re!, cdg, Q, beta
     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols
     type (tide_params) :: pars
     type(params)	:: P
     type(grid_dims):: GD

     integer ::      j, istat
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max
     integer, allocatable :: hp(:, :),up(:, :),vp(:, :)

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt
     integer :: iu(GD%nu), iv(GD%nv), ih(GD%nh)

     type (triplet) :: mat
     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr) :: mat_csr
     type (csr_cmplx) :: mat_csr_cmplx
     type (csc_cmplx) :: mat_csc_cmplx

	 type (csr) :: tmp_csr
	 type (csr_cmplx) :: tmp_csr_cmplx
	 type (triplet) :: tempmat
     integer ::      nmat, filled
	 type (triplet) :: DUU, DVU, DUV, DVV
     complex(cwp), allocatable :: Qu(:), Qv(:)!, tmp(:), Du(:), Dv(:)
	 real(wp), allocatable :: KU(:), KV(:)
	 real(wp), allocatable :: N_u(:), N_v(:), dHdph(:), dHdta(:)
     integer, allocatable, dimension(:) :: ufactor,vfactor

     type (csr_cmplx) :: speye
     integer, allocatable :: bcdiag(:)
     complex(cwp), allocatable :: rhs(:)

     complex(cwp), allocatable :: uvh(:), u(:), v(:), h(:)

!==========================================================================================
! SHORTCUTS
!==========================================================================================
latP = P%latP
lonP = P%lonP
g = P%g
re = P%re

nph = GD%nph
nth = GD%nta
nu = GD%nu
nv = GD%nv
nh = GD%nh
np = GD%np
!==========================================================================================
call tostring_set(rfmt='F12.1')
write(*, '("====================================")')
write(*, '("Welcome to baro_solver_linear")')
write(*, '("====================================")')


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

!%==================================================
!% Load the sparse matrix generated by baro_uvhT_mat
!%==================================================
     call load_alloc_sparse(mat, dir_mats // 'temp/mat_init.dat')
     call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')

! Convert mat from COO to CSR:
     call coo2csr(mat, mat_csr, P%lib)
     call dealloc_sparse(mat)

!  now transform into cmplx matrix
     call init_sparse_vals(mat_csr_cmplx,mat_csr%indi,mat_csr%indj,mat_csr%vals, &
                                mat_csr%ni, mat_csr%nj, mat_csr%nz)
      call dealloc_sparse(mat_csr)
!%==================================================

      call disp('Solving ' // tostring(mat_csr_cmplx%ni) //' x '// tostring(mat_csr_cmplx%nj)//' system with '// &
                 tostring(mat_csr_cmplx%nz)// ' entries for:' )
!==========================================================================================
!    %========================
!    % add internal tide drag ! should be incorporated into write_baro_mats
!    %========================
    if (P%itd_scheme == 1) then

		call load_alloc_vector(KU, dir_cols // 'KU.dat')
		call load_alloc_vector(KV, dir_cols // 'KV.dat')
		call load_alloc_sparse(DUU, dir_mats // 'DUU.dat')
		call load_alloc_sparse(DUV, dir_mats // 'DUV.dat')
		call load_alloc_sparse(DVU, dir_mats // 'DVU.dat')
		call load_alloc_sparse(DVV, dir_mats // 'DVV.dat')
		allocate(Qu(nu), Qv(nv), stat = istat)

		call calc_itd_prefactor(ufactor, vfactor, N_u, N_v, dHdph, dHdta, cpt, P, nu, nv, dir_cols, dir_grid)
		!***********************************************
		!	Add terms into the matrix
		!***********************************************
		filled = 0
		nmat = DUU%nz + DUV%nz + DVU%nz + DVV%nz
		call init_sparse_0(tempmat, np,np, nmat)

		! 1) DUU
		call left_mult_diag(DUU, ufactor*KU/pars%omega0, nu)
		call concatenate_sparse(tempmat, DUU, 0, 0, filled)
		! 2) DVU
		call left_mult_diag(DVU, ufactor*KU/pars%omega0, nu)
		call concatenate_sparse(tempmat, DVU, 0, nu, filled)
		! 3) DUV
		call left_mult_diag(DUV, vfactor*KV/pars%omega0, nv)
		call concatenate_sparse(tempmat, DUV, nu, 0, filled)
		! 4) DVV
		call left_mult_diag(DVV, vfactor*KV/pars%omega0, nv)
		call concatenate_sparse(tempmat, DVV, nu, nu, filled)

!		Convert the tempmat from COO to CSR:
	    call coo2csr(tempmat, tmp_csr, P%lib)
	    call init_sparse_vals(tmp_csr_cmplx,tmp_csr%indi,tmp_csr%indj,tmp_csr%vals, &
			                           tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
!	Save for testing
!		call save_vector( ufactor*KU/pars%omega0, dir_cols // 'tmp_u.dat')
!		call save_vector( vfactor*KV/pars%omega0, dir_cols // 'tmp_v.dat')
!     Save the final version of the matrix
!		call csr2coo(tmp_csr_cmplx, mat_coo_cmplx, P%lib)
!		call save_sparse(mat_coo_cmplx, dir_mats // 'mat1_' // cpt //'.dat')
!		call dealloc_sparse(mat_coo_cmplx)

	    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
	    call dealloc_sparse(tmp_csr)
	    call dealloc_sparse(tmp_csr_cmplx)

		call dealloc_sparse(DUU)
		call dealloc_sparse(DUV)
		call dealloc_sparse(DVU)
		call dealloc_sparse(DVV)
		deallocate(KU, KV, ufactor, vfactor)

	end if
!==========================================================================================

!     Save the final version of the matrix
!      call csr2coo(mat_csr_cmplx, mat_coo_cmplx, P%lib)
!      call save_sparse(mat_coo_cmplx, dir_mats // 'mat_' // cpt //'.dat')
!      call dealloc_sparse(mat_coo_cmplx)

!%=======================
!% calculate rhs forcing
!%=======================
     call baro_rhs(rhs, cpt, P, GD, dir_cols, dir_grid, dir_mats)
     call save_vector(rhs, dir_mats // 'rhs_' // cpt //'.dat')

!  %=============
!  % impose bcs:
!  %=============
write(*, '("Implementing boundary conditions...")', advance = 'no')
!mat=bcmat*mat; !fill the row of mat with zeroes for boundary variables
!rhs=bcmat*rhs; ! bcmat is already the way it should be!!!!!
     call dia_csr_mul(real(bcdiag,wp), mat_csr_cmplx, skit)
     rhs = rhs * bcdiag
     deallocate(bcdiag)
print *, "............. COMPLETE"

!==========================================================================================
!     mat = mat-i*omega0*speye(np)
      speye = init_speye(np, cmplx(0., -pars%omega0, kind=cwp) )
      call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), speye, mkl)

!  %=======
!  % solve
!  %=======

  write(*, '(" ", a, "... " )', advance='no') cpt
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

if (P%solver .eq. 'umfpack') then
!     USE UMFPACK TO SOLVE THE SYSTEM

      ! convert from csr to csc for umfpack solver (csc is convert to 0-based array)
        call csr2csc(mat_csr_cmplx, mat_csc_cmplx, P%lib)
        call dealloc_sparse(mat_csr_cmplx)
!                                           messages,filesave,loadsym,savenum
        call umfpack_unsym(mat_csc_cmplx%ni, mat_csc_cmplx, rhs, uvh, .false.,1,.false., .false.)
        call dealloc_sparse(mat_csc_cmplx)

elseif (P%solver .eq. 'pardiso') then
!     USE PARDISO TO SOLVE THE SYSTEM

	  call pardiso_unsym(mat_csr_cmplx%ni, mat_csr_cmplx, rhs, uvh, P, ccpt, .false., dir_mats // 'temp/') ! dir where factors can be saved
      call dealloc_sparse(mat_csr_cmplx)

else
      print *, "You must choose a valid solver: pardiso or umfpack"
      stop
end if

  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )
  write(*, '("done in ", f5.1, "s.")') T2-T1

!  %=================
!  % output solution
!  %=================
    allocate(u(nu), v(nv), h(nh), stat = istat)
    u = uvh(iu)
    v = uvh(iv)
    h = uvh(ih)
    deallocate(uvh)

     call save_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')

     call save_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
     call save_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

	deallocate(u, v, h)

enddo

!=========================================================================================
!=================================DONE====================================================
!=========================================================================================
write(*, '("========================================================================")')
write(*, '("========================================================================")')

do ccpt = 1, ncpts

     cpt=cpts(2*ccpt-1:2*ccpt)

     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
     call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

     write(*,*) cpt, " max tide height: ", maxval(abs(h))

!  if flags.graphics == 1
!
!    %==============================
!    % plot global tides and errors
!    %==============================

     if (ccpt == 1) then

         call load_alloc_matrix(hp, dir_cols // 'hp.dat')
         call load_alloc_matrix(up, dir_cols // 'up.dat')
         call load_alloc_matrix(vp, dir_cols // 'vp.dat')

		call save_sol4plot_cmplx(nph, nth, nh, hp, h, dir_sols // 'sol_h_'//cpt//'.dat')
		call save_sol4plot_cmplx(nph, nth, nu, up, u, dir_sols // 'sol_u_'//cpt//'.dat')
		call save_sol4plot_cmplx(nph, nth, nv, vp, v, dir_sols // 'sol_v_'//cpt//'.dat')

		deallocate(hp, up, vp)

     endif

	deallocate(u,v,h)

enddo

end subroutine baro_solver_linear
!==========================================================================================

!==========================================================================================
subroutine baro_rhs(rhs, cpt, P, GD, dir_cols, dir_grid, dir_mats)

    implicit none
!
!% calculates the complex barotropic tidal forcing
!% for the given tidal component.

     character(len=2), intent(in)   :: cpt
     type(params), intent(in)		:: P
     type(grid_dims), intent(in)	:: GD

     integer  :: nu, nv, nh, coor
     real(wp) :: latP, lonP

     type (triplet) :: GHU, GHV

     integer ::      j
     integer ::      istatus
     integer, allocatable :: 	 iu(:), iv(:)
!     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     complex(cwp), allocatable, dimension(:) :: rhs, heq
!     real(wp), allocatable, dimension(:) :: tmp2, ph_h, ph_u

     character(len = *) :: dir_grid, dir_cols, dir_mats

coor = P%coor
latP =  P%latP
lonP = P%lonP

nu = GD%nu
nv = GD%nv
nh = GD%nh


	allocate(iu(nu), iv(nv), stat=istatus)
     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (j), j=1,nv ) /)
!%==============================================================
!% Load the necessary matrices:
!% GHU, GHV
!%==============================================================
     call load_alloc_sparse(GHU, dir_mats // 'GHU.dat')
     call load_alloc_sparse(GHV, dir_mats // 'GHV.dat')
!*************************************************************************
     allocate(rhs(nu+nv+nh), stat = istatus)
     rhs = 0.

     call calc_heq_h(heq, cpt, nh, latP, lonP, coor, dir_cols, dir_grid)
!     call save_vector(heq, dir_cols // 'heq.dat')

     call coo_vec_mul(GHU, heq, P%lib, rhs(1:nu))
     call dealloc_sparse(GHU)

     call coo_vec_mul(GHV, heq, P%lib, rhs(nu+1:nu+nv))
     call dealloc_sparse(GHV)

!     call save_vector(rhs, dir_cols // 'rhs.dat')

end subroutine baro_rhs
!==========================================================================================

!==========================================================================================
subroutine prepare_init_guess(P, GD, out_dir, P_file, GD_file, tidal_dir, dir_grid, dir_cols, dir_mats, dir_sols)

implicit none

     type(params)	:: P
     type(grid_dims), intent(in):: GD
     character(len = *)		:: out_dir,tidal_dir, dir_sols, dir_grid, dir_cols, dir_mats
     character(len = *)		:: P_file, GD_file

     type(params)			:: P_load
     type(grid_dims)		:: GD_load

	 logical			:: dir_e
	 logical, allocatable	:: tpxo_e(:),  prev_sol_e(:)
     character(len=2)	:: cpt
     integer            :: j, ccpt, ncpts, istat
     character(len=100)	:: sal_cpts

     complex(cwp), allocatable	:: u(:), v(:), h(:), Du_bbl(:), Dv_bbl(:), hsal(:), hsal_hg_sol(:,:)
     complex(cwp), allocatable	:: u_sol(:), v_sol(:), h_sol(:), Du_bbl_sol(:), Dv_bbl_sol(:), hsal_sol(:)
     complex(cwp), allocatable	:: h_hg_sol(:,:), u_ug_sol(:,:), v_vg_sol(:,:), Du_ug_bbl_sol(:,:), Dv_vg_bbl_sol(:,:)
     complex(cwp), allocatable	:: h_hg(:,:), u_ug(:,:), v_vg(:,:), Du_ug_bbl(:,:), Dv_vg_bbl(:,:), hsal_hg(:,:)
     integer, allocatable		:: hp(:,:), up(:,:), vp(:,:)
     integer, allocatable		:: hp_sol(:,:), up_sol(:,:), vp_sol(:,:)
     real(wp), allocatable		:: beta(:), beta_sol(:), beta_hg(:,:), beta_hg_sol(:,:), beta0(:)
     real(wp), allocatable		:: ph_vg(:), th_ug(:), ph_vg_sol(:), th_ug_sol(:)
     real(wp), allocatable		:: ph_ug(:), th_vg(:), ph_ug_sol(:), th_vg_sol(:)
     real(wp), allocatable		:: H_h(:), H_u(:), H_v(:)

     complex(cwp), allocatable	:: hg_tpxo(:,:), ug_tpxo(:,:), vg_tpxo(:,:)
     real(wp), allocatable		:: lon(:), lat(:)
     integer					:: nlon, nlat

!     character(len = *)		:: dir_grid, dir_cols, dir_mats, N_data_dir

write(*, '("----------------------------------------------------------------------------------------")')

if (P%load_sol == 1) then

	write(*, '("Prepare the initial guess baro solution: use ", a)') trim(P%gridid_sol)

	! check that the directory exists
		inquire( file=out_dir//trim(P%gridid_sol)//'/.', exist=dir_e )
		if  (.not. dir_e ) then
			write(*,'(a)') 'WARNING: Directory /'// out_dir//trim(P%gridid)//' for the previously generated solution doesn''t exist.'
			write(*,*) "Solving with zero initial guess"
			write(*,*) ""
			P%load_sol = 0
			return
!			stop
		end if
	! The directory with the previously generated solution exists...

	!	Load grid and parameters file
	call control_file(out_dir // trim(P%gridid_sol) // '/' // P_file, P_load)
	!	Ad-hoc fix: SAL directory in "/sols/temp/' is normaly cleaned up. So, don't do interpolation
		if (P_load%sal_scheme >=3) then
			P_load%sal_scheme = 0
		endif

	call read_GD(out_dir // trim(P%gridid_sol) // '/' // 'global/grid/'//GD_file, GD_load)

	ncpts=len(trim(P%cpts))/2
	! Index cpts that were computed in the "previous run"
		allocate(prev_sol_e(ncpts), stat=istat)
			prev_sol_e=.false.
			do ccpt = 1, ncpts
		    	cpt=P%cpts(2*ccpt-1:2*ccpt)
				prev_sol_e(ccpt) = INDEX(P_load%cpts, cpt) > 0
			end do
			if  (.not. (count(prev_sol_e)>0) ) then
				write(*,'(a)') 'WARNING: Directory /'// out_dir//trim(P%gridid)//' does not contain required cpts.'
				write(*,*) "Solving with zero initial guess"
				write(*,*) ""
				P%load_sol = 0
				return
			endif

			! make a list of cpts in both P_load and P
			j=0
			do ccpt = 1, ncpts
				if (prev_sol_e(ccpt)) then
					j = j+1
					sal_cpts(2*j-1:2*j) = P%cpts(2*ccpt-1:2*ccpt)
		    	endif
			end do


	! check if grids are identical. If yes, just link the corresponding files over to the new directory
	if (cmpr_GD(GD_load, GD)) then
		!% Copy previously generated files
		write(*, '("Identical grids. Present solution is copied and taken as initial guess.")')
	     call system('cp ' // out_dir//trim(P%gridid_sol)//'/'//'global/sols/'//'*m0.dat ' // dir_sols) ! u, v, h -solutions
	     call system('cp ' // out_dir//trim(P%gridid_sol)//'/'//'global/sols/'//'*m0_drag.dat ' // dir_sols) ! corresponding BBL drag
	     call system('cp ' // out_dir//trim(P%gridid_sol)//'/'//'global/sols/temp/'//'*sal.dat ' // dir_sols//'/temp/') ! corresponding SAL terms
	    ! if SAL terms were not computed or was cleaned in the previous run, then compute them now
	    if ((P%sal_scheme >= 2).and.(P_load%sal_scheme < 3) .and.(P%itd_scheme /= 20) .or. (P_load%cleanup>=1)) then
	    	! only run calc_sal for cpts present in both P_load and P: written in sal_cpts
	    	j = count(prev_sol_e)
			allocate(beta0(j), stat = istat)
			beta0 = P%beta0
	    	call calc_sal(j, sal_cpts(1:2*j), P, GD, dir_grid, dir_cols, dir_sols, beta0)
	    endif
	else ! do bilinear interpolation. WARNING: only works if coordinates of the rotated poles are the same.
		! check this
		if  (.not. ( (GD_load%lonP == GD%lonP).and.(GD_load%latP == GD%latP) ) ) then
			write(*,'(a)') 'WARNING: using a solution on a grid with a different rotated pole is NOT implemented.'
			write(*,'(a)') "Solving with zero initial guess"
			P%load_sol = 0
			return
		end if
		!    %======================================
		!    % load grids: interpolating from and to
		!    %======================================
!			write(*, '("1) loading grids for interpolation: ")')
			!***************************** Solution grid *****************************
			! load hp, up, vp
		     call load_alloc_matrix(up_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/cols/' // 'up.dat')
		     call load_alloc_matrix(vp_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/cols/' // 'vp.dat')
		     call load_alloc_matrix(hp_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/cols/' // 'hp.dat')
		     ! grid coords ph_vg, ph_vg, th_ug, th_vg
		     call load_alloc_vector(ph_vg_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/grid/' // 'ph_vg.dat')
		     call load_alloc_vector(th_ug_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/grid/' // 'th_ug.dat')
		     call load_alloc_vector(ph_ug_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/grid/' // 'ph_ug.dat')
		     call load_alloc_vector(th_vg_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/grid/' // 'th_vg.dat')

		    allocate(u_ug_sol(GD_load%nph, GD_load%nta), v_vg_sol(GD_load%nph, GD_load%nta), &
		    		 h_hg_sol(GD_load%nph, GD_load%nta), stat = istat)
		    u_ug_sol = 0
		    v_vg_sol = 0
		    h_hg_sol = 0
			if (P_load%sal_scheme >= 3) then
				allocate(hsal_hg_sol(GD_load%nph, GD_load%nta), beta_hg_sol(GD_load%nph, GD_load%nta),  stat = istat)
			    hsal_hg_sol = 0
			    beta_hg_sol = 0
			endif
			if (P_load%fr_scheme >= 3) then
				allocate(Du_ug_bbl_sol(GD_load%nph, GD_load%nta), Dv_vg_bbl_sol(GD_load%nph, GD_load%nta),  stat = istat)
			    Du_ug_bbl_sol = 0
			    Dv_vg_bbl_sol = 0
			endif

		    !***************************** Interpolation grid *****************************
			! load hp, up, vp
		     call load_alloc_matrix(up, dir_cols // 'up.dat')
		     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
		     call load_alloc_matrix(hp, dir_cols // 'hp.dat')
		    ! grid coords ph_vg, ta_h, H_hg
		    call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
		    call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
		    call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
		    call load_alloc_vector(th_vg, dir_grid // 'th_vg.dat')
			call load_alloc_vector(H_h, dir_cols // 'H_h.dat')
			call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
			call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
			! allocate interpolated u and v and h
			allocate(u(GD%nu),v(GD%nv), h(GD%nh), stat = istat)
			if (P%sal_scheme >= 3) then
				allocate(hsal(GD%nh), beta(GD%nh),  stat = istat)
			endif
			if (P%fr_scheme >= 3) then
				allocate(Du_bbl(GD%nu), Dv_bbl(GD%nv),  stat = istat)
			endif


		if (P%sal_scheme >= 2) then
			allocate(beta0(ncpts), stat = istat)
			beta0 = P%beta0
		endif

		do ccpt = 1, ncpts
		  if (prev_sol_e(ccpt)) then
	    	cpt=P%cpts(2*ccpt-1:2*ccpt)

			write(*, '(i2, a)', advance = 'no') ccpt, ") "//cpt//": a) load u, v, h"
			! Load u, v, h -solutions
			call load_alloc_vector(u_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/' // cpt // '_u' // '_m0' // '.dat')
			call load_alloc_vector(v_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/' // cpt // '_v' // '_m0' // '.dat')
			call load_alloc_vector(h_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/' // cpt // '_h' // '_m0' // '.dat')
			! Write solutions on a grid for interpolation
			do j = 1,GD_load%nh
				h_hg_sol(hp_sol(j,1),hp_sol(j,2)) = h_sol(j)
			enddo
			do j = 1,GD_load%nu
				u_ug_sol(up_sol(j,1),up_sol(j,2)) = u_sol(j)
			enddo
			do j = 1,GD_load%nv
				v_vg_sol(vp_sol(j,1),vp_sol(j,2)) = v_sol(j)
			enddo

		  if ((P%sal_scheme >= 3).and.(P_load%sal_scheme >= 3) .and.(P%itd_scheme /= 20)) then
			write(*, '(a)', advance = 'no') ", betasal, hsal"
			! Load Du, Dv for BBL
			call load_alloc_vector(beta_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/temp/' // cpt // '_betasal' // '.dat')
			call load_alloc_vector(hsal_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/temp/' // cpt // '_hsal' // '.dat')
			! Write solutions on a grid for interpolation
			do j = 1,GD_load%nh
				beta_hg_sol(hp_sol(j,1),hp_sol(j,2)) = beta_sol(j)
				hsal_hg_sol(hp_sol(j,1),hp_sol(j,2)) = hsal_sol(j)
			enddo
		  endif

		  if ((P%fr_scheme >= 3).and.(P_load%fr_scheme >= 3)) then
			write(*, '(a)', advance = 'no') ", Du_m0_drag, Dv_m0_drag on a grid"
			! Load Du, Dv for BBL
			call load_alloc_vector(Du_bbl_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/' // cpt // '_Du_m0_drag' // '.dat')
			call load_alloc_vector(Dv_bbl_sol, out_dir // trim(P%gridid_sol) // '/' // 'global/sols/' // cpt // '_Dv_m0_drag' // '.dat')
			! Write solutions on a grid for interpolation
			do j = 1,GD_load%nu
				Du_ug_bbl_sol(up_sol(j,1),up_sol(j,2)) = Du_bbl_sol(j)
			enddo
			do j = 1,GD_load%nv
				Dv_vg_bbl_sol(vp_sol(j,1),vp_sol(j,2)) = Dv_bbl_sol(j)
			enddo
		  endif

			! interpolate onto the new grid
			write(*, '(a)', advance = 'no') "; b) interpolating: "
			write(*, '(a)', advance = 'no') "u, "
			call add_border_for_interp(GD_load%nph, GD_load%nta,u_ug_sol, abs(u_ug_sol)==0)	! for better interpolation along the land border
	    	call bilinear_2d(GD_load%nph, GD_load%nta, ph_ug_sol, th_ug_sol, u_ug_sol, GD%nph, GD%nta, ph_ug, th_ug, u_ug, P%lib)
	    	write(*, '(a)', advance = 'no') "v, "
	    	call add_border_for_interp(GD_load%nph, GD_load%nta, v_vg_sol, abs(v_vg_sol)==0)	! for better interpolation along the land border
	    	call bilinear_2d(GD_load%nph, GD_load%nta, ph_vg_sol, th_vg_sol, v_vg_sol, GD%nph, GD%nta, ph_vg, th_vg, v_vg, P%lib)
	    	write(*, '(a)', advance = 'no') "h"
	    	call add_border_for_interp(GD_load%nph, GD_load%nta, h_hg_sol, abs(h_hg_sol)==0)	! for better interpolation along the land border
	    	call bilinear_2d(GD_load%nph, GD_load%nta, ph_vg_sol, th_ug_sol, h_hg_sol, GD%nph, GD%nta, ph_vg, th_ug, h_hg, P%lib)

		if ((P%sal_scheme >= 3).and.(P_load%sal_scheme >= 3) .and.(P%itd_scheme /= 20)) then
			write(*, '(a)', advance = 'no') ", beta, "
		call add_border_for_interp(GD_load%nph, GD_load%nta, beta_hg_sol, abs(beta_hg_sol)==0)	! for better interpolation along the land border
	    call bilinear_2d(GD_load%nph, GD_load%nta, ph_vg_sol, th_ug_sol, beta_hg_sol, GD%nph, GD%nta, ph_vg, th_ug, beta_hg, P%lib)
			write(*, '(a)', advance = 'no') "hsal"
	    call add_border_for_interp(GD_load%nph, GD_load%nta, hsal_hg_sol, abs(hsal_hg_sol)==0)	! for better interpolation along the land border
	    call bilinear_2d(GD_load%nph, GD_load%nta, ph_vg_sol, th_ug_sol, hsal_hg_sol, GD%nph, GD%nta, ph_vg, th_ug, hsal_hg, P%lib)
	    endif

		if ((P%fr_scheme >= 3).and.(P_load%fr_scheme >= 3)) then
			write(*, '(a)', advance = 'no') ", Du_m0_drag, "
		call add_border_for_interp(GD_load%nph, GD_load%nta, Du_ug_bbl_sol, abs(Du_ug_bbl_sol)==0)	! for better interpolation along the land border
	    call bilinear_2d(GD_load%nph, GD_load%nta, ph_ug_sol, th_ug_sol, Du_ug_bbl_sol, GD%nph, GD%nta, ph_ug, th_ug, Du_ug_bbl, P%lib)
			write(*, '(a)', advance = 'no') "Dv_m0_drag"
		call add_border_for_interp(GD_load%nph, GD_load%nta, Dv_vg_bbl_sol, abs(Dv_vg_bbl_sol)==0)	! for better interpolation along the land border
	    call bilinear_2d(GD_load%nph, GD_load%nta, ph_vg_sol, th_vg_sol, Dv_vg_bbl_sol, GD%nph, GD%nta, ph_vg, th_vg, Dv_vg_bbl, P%lib)
	    endif
		!call save_matrix(h_hg, dir_grid // 'N_hg.dat')

		 ! write u, v, h in columns
		    do j = 1,GD%nu
		    	u(j) = u_ug(up(j,1),up(j,2))
		    enddo
		    where (.not.(H_u > 0)) u = 0

		    do j = 1,GD%nv
		    	v(j) = v_vg(vp(j,1),vp(j,2))
		    enddo
		    where (.not.(H_v > 0)) v = 0

			if ((P%fr_scheme >= 3).and.(P_load%fr_scheme >= 3) .and.(P%itd_scheme /= 20)) then
		    	do j = 1,GD%nu
	    			Du_bbl(j) = Du_ug_bbl(up(j,1),up(j,2))
			    enddo
			    do j = 1,GD%nv
		    		Dv_bbl(j) = Dv_vg_bbl(vp(j,1),vp(j,2))
		    	enddo

		    	where (.not.(H_u > 0)) Du_bbl = 0
		    	where (.not.(H_v > 0)) Dv_bbl = 0
		    endif

		    do j = 1,GD%nh
 				h(j) = h_hg(hp(j,1),hp(j,2))
 			enddo
 			where (.not.(H_h > 0)) h = 0

		    if ((P%sal_scheme >= 3).and.(P_load%sal_scheme >= 3)) then
		    	do j = 1,GD%nh
	    		hsal(j) = hsal_hg(hp(j,1),hp(j,2))
	    		beta(j) = beta_hg(hp(j,1),hp(j,2))
			    enddo

		    	where (.not.(H_h > 0)) hsal = 0
		    	where (.not.(H_h > 0)) beta = 0
		    endif

			deallocate(u_ug, v_vg, h_hg)

			! Save u, v, h -solutions
		     call save_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
		     call save_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
		     call save_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

			if ((P%sal_scheme >= 3).and.(P_load%sal_scheme >= 3) .and.(P%itd_scheme /= 20)) then
				deallocate(hsal_hg, beta_hg)
				call save_vector(hsal, dir_sols//'temp/' // cpt // '_hsal' // '.dat')
		    	call save_vector(beta, dir_sols//'temp/' // cpt // '_beta_sal' // '.dat')
			endif

			if ((P%fr_scheme >= 3).and.(P_load%fr_scheme >= 3)) then
				deallocate(Du_ug_bbl, Dv_vg_bbl)
				call save_vector(Du_bbl, dir_sols // cpt // '_Du_m0_drag' // '.dat')
		    	call save_vector(Dv_bbl, dir_sols // cpt // '_Dv_m0_drag' // '.dat')
			endif

			write(*, *) "... DONE."
		  endif
	    enddo

		! if the previous run does not contain all the
			do ccpt = 1, ncpts
	    		cpt=P%cpts(2*ccpt-1:2*ccpt)

				if (.not. prev_sol_e(ccpt)) then
					write(*, '(i2, a)', advance = 'no') ccpt, ") "//cpt//": Previously generated solution is NOT available. Set u, v, h = 0"
		  			u=0; v=0; h=0;
					! Save u, v, h -solutions
				     call save_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
				     call save_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
				     call save_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

					if ((P%sal_scheme >= 3)) then
						hsal = 0; beta = P%beta0
						call save_vector(hsal, dir_sols//'temp/' // cpt // '_hsal' // '.dat')
				    	call save_vector(beta, dir_sols//'temp/' // cpt // '_betasal' // '.dat')
					endif

					if ((P%fr_scheme >= 3).and.(P_load%fr_scheme >= 3)) then
						Du_bbl = 0; Dv_bbl=0
						call save_vector(Du_bbl, dir_sols // cpt // '_Du_m0_drag' // '.dat')
				    	call save_vector(Dv_bbl, dir_sols // cpt // '_Dv_m0_drag' // '.dat')
					endif

					write(*, *) "... DONE."

				endif
			end do

		! if SAL terms were not computed in the previous run, then compute them now
	    if ((P%sal_scheme >= 2).and.(P_load%sal_scheme < 3) .and.(P%itd_scheme /= 20)) then
	    	! only run calc_sal for cpts present in both P_load and P: written in sal_cpts
	    	j = count(prev_sol_e)
	    	call calc_sal(j, sal_cpts(1:2*j), P, GD, dir_grid, dir_cols, dir_sols, beta0)
	    endif
		! if BBL friction terms were not computed in the previous run, then compute them now
	    if ((P%fr_scheme >= 3).and.(P_load%fr_scheme < 3)) then
	    	call calc_dhat(trim(P%cpts), GD%nu, GD%nv, P, dir_cols, dir_mats, dir_sols)
	    endif
	!==========================================================================================

		deallocate(u, v,h)
		deallocate(u_ug_sol, v_vg_sol, h_hg_sol)

!    deallocate(hp, ph_vg, th_ug, N_h)
	endif
elseif (P%load_sol == 2) then

	call tostring_set_factory()

	write(*, '("Prepare the initial guess baro solution: use ", a)') trim(P%tpxo_dir)

	! check that the directory exists
		inquire( file=tidal_dir//trim(P%tpxo_dir)//'/.', exist=dir_e )
		if  (.not. dir_e ) then
			write(*,'(a)') 'WARNING: Directory /'// tidal_dir//trim(P%tpxo_dir)//' for the TPXO solution doesn''t exist.'
			write(*,*) "Solving with zero initial guess"
			write(*,*) ""
			P%load_sol = 0
			return
		else
			ncpts=len(trim(P%cpts))/2
			allocate(tpxo_e(ncpts), stat=istat)
			tpxo_e=.false.
			do ccpt = 1, ncpts
		    	cpt=P%cpts(2*ccpt-1:2*ccpt)
				inquire( file=tidal_dir//trim(P%tpxo_dir)//'/h_'//cpt//'_rot_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat', &
																												exist=tpxo_e(ccpt) )
			end do
			if  (.not. (count(tpxo_e)>0) ) then
				write(*,'(a)') 'WARNING: Files '//trim(P%tpxo_dir)//'/h_**_rot_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat'// &
										' for the TPXO solution do not exist.'
				write(*,*) "Solving with zero initial guess"
				write(*,*) ""
				P%load_sol = 0
				return
			endif
		endif

			! make a list of cpts in both P and TPXO
			j=0
			do ccpt = 1, ncpts
				if (tpxo_e(ccpt)) then
					j = j+1
					sal_cpts(2*j-1:2*j) = P%cpts(2*ccpt-1:2*ccpt)
		    	endif
			end do

	! The TPXO directory and files do exist...
		!    %======================================
		!    % load grids: interpolating from and to
		!    %======================================
		    !***************************** Interpolation grid *****************************
			! load hp, up, vp
		     call load_alloc_matrix(up, dir_cols // 'up.dat')
		     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
		     call load_alloc_matrix(hp, dir_cols // 'hp.dat')
		    ! grid coords ph_vg, ta_h, H_hg
		    call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
		    call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
		    call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
		    call load_alloc_vector(th_vg, dir_grid // 'th_vg.dat')
			call load_alloc_vector(H_h, dir_cols // 'H_h.dat')
			call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
			call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
			! allocate interpolated u and v and h
			allocate(u(GD%nu),v(GD%nv), h(GD%nh), stat = istat)
			if (P%sal_scheme >= 3) then
				allocate(hsal(GD%nh), beta(GD%nh),  stat = istat)
			endif
			if (P%fr_scheme >= 3) then
				allocate(Du_bbl(GD%nu), Dv_bbl(GD%nv),  stat = istat)
			endif


		ncpts=len(trim(P%cpts))/2

		if (P%sal_scheme >= 2) then
			allocate(beta0(ncpts), stat = istat)
			beta0 = P%beta0
		endif

		do ccpt = 1, ncpts
	    	cpt=P%cpts(2*ccpt-1:2*ccpt)

		  if (tpxo_e(ccpt)) then
			write(*, '(i2, a)', advance = 'no') ccpt, ") "//cpt//": a) load u, v, h"
			!***************************** TPXO grid *****************************
			! Load TPXO u, v, h -solutions on a grid
			call load_alloc_topo(tidal_dir//trim(P%tpxo_dir)//'/h_'//cpt//'_rot_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat', &
								nlon, nlat, lon, lat, hg_tpxo, .false.)
			call load_alloc_topo(tidal_dir//trim(P%tpxo_dir)//'/u_'//cpt//'_rot_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat', &
								nlon, nlat, lon, lat, ug_tpxo, .false.)
			call load_alloc_topo(tidal_dir//trim(P%tpxo_dir)//'/v_'//cpt//'_rot_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat', &
								nlon, nlat, lon, lat, vg_tpxo, .false.)

			! interpolate onto the new grid
			write(*, '(a)', advance = 'no') "; b) interpolating: "
			write(*, '(a)', advance = 'no') "u, "
			call add_border_for_interp(nlon, nlat, ug_tpxo, abs(hg_tpxo)==0)	! for better interpolation along the land border
	    	call bilinear_2d(nlon, nlat, lon, lat, ug_tpxo, GD%nph, GD%nta, ph_ug, th_ug, u_ug, P%lib)
	    	write(*, '(a)', advance = 'no') "v, "
	    	call add_border_for_interp(nlon, nlat, vg_tpxo, abs(hg_tpxo)==0)	! for better interpolation along the land border
	    	call bilinear_2d(nlon, nlat, lon, lat, vg_tpxo, GD%nph, GD%nta, ph_vg, th_vg, v_vg, P%lib)
	    	write(*, '(a)', advance = 'no') "h"
	    	call add_border_for_interp(nlon, nlat, hg_tpxo, abs(hg_tpxo)==0)	! for better interpolation along the land border
	    	call bilinear_2d(nlon, nlat, lon, lat, hg_tpxo, GD%nph, GD%nta, ph_vg, th_ug, h_hg, P%lib)

		 ! write u, v, h in columns
		    do j = 1,GD%nu
		    	u(j) = u_ug(up(j,1),up(j,2))
		    enddo
		    where (.not.(H_u > 0)) u = 0

		    do j = 1,GD%nv
		    	v(j) = v_vg(vp(j,1),vp(j,2))
		    enddo
		    where (.not.(H_v > 0)) v = 0

		    do j = 1,GD%nh
 				h(j) = h_hg(hp(j,1),hp(j,2))
 			enddo
 			where (.not.(H_h > 0)) h = 0

!			deallocate(u_ug, v_vg, h_hg)
		  else
			write(*, '(i2, a)', advance = 'no') ccpt, ") "//cpt//": TPXO solution is NOT available. Set u, v, h = 0"
		  	u=0; v=0; h=0
		  endif
			! Save u, v, h -solutions
		     call save_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
		     call save_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
		     call save_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

			write(*, *) "... DONE."
	    enddo

		! compute SAL terms now
	    if ((P%sal_scheme >= 2) .and.(P%itd_scheme /= 20)) then
	    	! only run calc_sal for cpts present in both P_load and P: written in sal_cpts
	    	j = count(tpxo_e)
	    	call calc_sal(j, sal_cpts(1:2*j), P, GD, dir_grid, dir_cols, dir_sols, beta0)

	    	! ad-hoc correction. Beta can be NAN at Black sea (and such). Correct
			do ccpt = 1, ncpts
		    	cpt=P%cpts(2*ccpt-1:2*ccpt)

				call load_vector(beta, dir_sols //'temp/'// cpt // '_betasal' // '.dat')
		    	where (isnan(beta)) beta = P%beta0
				call save_vector(beta, dir_sols //'temp/'// cpt // '_betasal' // '.dat')
		    enddo
	    endif

		! compute BBL friction terms now
	    if (P%fr_scheme >= 3) then
	    	call calc_dhat(trim(P%cpts), GD%nu, GD%nv, P, dir_cols, dir_mats, dir_sols)
	    endif
	!==========================================================================================

		deallocate(u, v,h)
		deallocate(ug_tpxo, vg_tpxo, hg_tpxo)

endif

write(*, '("----------------------------------------------------------------------------------------")')

end subroutine prepare_init_guess

!**********************************************************************************************************
!**********************************************************************************************************

end module baro_solver_mod



