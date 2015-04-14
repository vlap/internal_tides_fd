module generate_matrices

     use precisions, only: dp, wp, cwp
     use dispmodule
     use my_trigs
     use my_sparse
     use my_sparse_aggregate
     use save_load
     use generate_grid
     use itd
	 use control

     contains

!    %========================================================
!    % Generate aggregate matrices
!    %========================================================
!==========================================================================================
subroutine write_baro_mats(N_data_dir, dir_cols, dir_grid, dir_mats, cpts, P, GD)

     character(len=*) :: dir_cols, dir_grid, dir_mats
     character(len = *) :: N_data_dir
     type(params), intent(in)	 :: P
     type(grid_dims), intent(in) :: GD
     real(wp)	  :: cdg, udg, Q
     integer :: nu, nv, nh, np, cooruv, nph, nth
     real(wp)	:: re

     integer, allocatable :: up(:, :), vp(:, :)
     integer, allocatable :: iu(:), iv(:)!, ih(nh), ip(np)
     real(wp), allocatable ::  ta_u(:), ta_v(:), ph_vg(:), ta_ug(:)!, ta_h(:)

!     real    ::      T1, T2
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     integer ::      j, cu, cv
     integer ::      istat, nf!, statusj, statusval

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt
	 type (tide_params) :: pars

     real(wp), allocatable :: Qu(:), Qv(:)
     real(wp), allocatable :: H_u(:), H_v(:), H_h(:), KU(:), KV(:)
     real(wp), allocatable :: dHdph(:), dHdta(:), N_u(:), N_v(:)
	 real(wp), allocatable :: dHdph_prmtz(:), dHdta_prmtz(:), H_sm_sc_h(:) !, gradH_u_prmtz(:), gradH_v_prmtz(:)
	 real(wp), allocatable :: norm_u(:), norm_v(:), ratio_u(:), ratio_v(:)
     real(wp), allocatable :: norm_ug(:,:), norm_vg(:,:), ratio_ug(:,:), ratio_vg(:,:)
	 integer, allocatable    	:: mask_ug(:,:), mask_vg(:,:)
     real(wp), allocatable		:: metrics_ug(:, :), metrics_vg(:, :)

	 type (triplet)			:: tempmat
     integer				:: nmat, filled ! points how many elements out of nmat are already filled
	 integer, allocatable	:: ufactor(:), vfactor(:)
	 type (csr) :: tmp_csr

     type (triplet) :: GHU, GHV, JUH, JVH, DUU, DVU, DUV, DVV
     type (triplet) :: u2vf, v2uf, u2v, v2u, h2uddph, h2vddta

     type (triplet) :: mat
	 type (csr) :: mat_csr, mat_csr_cpt
     integer, pointer :: bcdiag(:)

!%==============================
write(*, '("Preparing matrices: ")', advance = 'no')

!%=============================================
! shortcuts
!%=============================================
!g = P%g
nph = GD%nph
nth = GD%nta

nu = GD%nu
nv = GD%nv
nh = GD%nh
re = P%re
np = nu + nv + nh
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
!%=============================================
!    %=========================================
!    % Martices for baro tides and not just ITs
!    %=========================================
	allocate(iu(nu),iv(nv), stat=istat)
     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (j), j=1,nv ) /)
!     ih = (/ ( (j), j=1,nh ) /)
!     ip = (/ ( (j), j=1,np ) /)

!    %==========================
!    % set up the basic matrices
!    %==========================
      call baro_uvhT_mat(nu, nv, nh, P, nmat, dir_cols, dir_mats)

!    %=========================================
!    % Check if baro solver will be incurred or
!	 % forced lin intern tides are calculated
!    %=========================================
if (P%itd_scheme /= 20) then
!%=============================================
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
!	Target:
!	nmat = v2uf%nz + u2vf%nz + h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz + nu + nv
!	After calling routine baro_uvhT_mat:
!	nmat = h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz+ nu + nv
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
	call load_alloc_sparse(v2uf, dir_mats // 'v2uf.dat')
	call load_alloc_sparse(u2vf, dir_mats // 'u2vf.dat')

	nmat = nmat + v2uf%nz + u2vf%nz
	filled = 0
	call init_sparse_0(mat, np,np, nmat)

!%===============
!	Add rotation
!%===============
!	1) v2uf
	v2uf%vals = - v2uf%vals ! multiply by -1
	call concatenate_sparse(mat, v2uf, 0, nu, filled)
	call dealloc_sparse(v2uf)

	write(*, '("v2uf, ")', advance = 'no')
!	2) u2vf
    call concatenate_sparse(mat, u2vf, nu, 0, filled)
	call dealloc_sparse(u2vf)

	write(*, '("u2vf, ")', advance = 'no')

!%==============================
!	Add already generated terms
!%==============================
!	Gravity (linear SAL is included if sal_scheme==1):
	call load_alloc_sparse(GHU, dir_mats // 'GHU.dat')
	call load_alloc_sparse(GHV, dir_mats // 'GHV.dat')

	if (P%sal_scheme == 1) then
		GHU%vals = (1 - P%beta0)*GHU%vals
		GHV%vals = (1 - P%beta0)*GHV%vals
	end if

	call concatenate_sparse(mat, GHU, 0, nu+nv, filled)
	call concatenate_sparse(mat, GHV, nu, nu+nv, filled)

	call dealloc_sparse(GHU)
	call dealloc_sparse(GHV)

!%===================================
!	Advection gradients (c.o. mass)
!%===================================
	call load_alloc_sparse(JUH, dir_mats // 'JUH.dat')
	call load_alloc_sparse(JVH, dir_mats // 'JVH.dat')

	call concatenate_sparse(mat, JUH, nu+nv, 0, filled)
	call concatenate_sparse(mat, JVH, nu+nv, nu, filled)

	call dealloc_sparse(JUH)
	call dealloc_sparse(JVH)

!%===================================
!	Linear bottom friction
!%===================================
	call load_alloc_vector(KU, dir_cols // 'KU.dat')
	call load_alloc_vector(KV, dir_cols // 'KV.dat')

	call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	call load_alloc_vector(H_v, dir_cols // 'H_v.dat')

!	TRANSPORT FORMULATION
	allocate(Qu(nu), Qv(nv), stat = istat)
		cdg = P%cd
	if (P%fr_scheme == 1) then
		udg = P%ubar
		Qu = cdg*udg
		Qv = cdg*udg
	elseif (P%fr_scheme == 2) then
		Q = P%Qbar
		Qu = cdg*Q/H_u
		Qv = cdg*Q/H_v
	else
		Qu = 0
		Qv = 0
	end if

	if (P%coor == 1) then
		call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
		call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')

		Qu=Qu * cosh(ta_u)
		Qv=Qv * cosh(ta_v)

		deallocate(ta_u, ta_v)
	endif
!      if (uvform == 2) then
        Qu=Qu/H_u
        Qv=Qv/H_v
!      endif


	call init_sparse(tempmat, iu, iu, KU*Qu, nu, nu, nu)
	call concatenate_sparse(mat, tempmat, 0, 0, filled)
	call dealloc_sparse(tempmat)

	call init_sparse(tempmat, iv, iv, KV*Qv, nv,nv, nv)
	call concatenate_sparse(mat, tempmat, nu, nu, filled)
	call dealloc_sparse(tempmat)
	write(*, '("Linear BBL, ")', advance = 'no')

	deallocate(Qu, Qv, H_u, H_v, KU, KV)

! Convert mat from COO to CSR:
     call coo2csr(mat, mat_csr, P%lib)
     call save_sparse(mat, dir_mats // 'temp/mat_init.dat')
     call dealloc_sparse(mat)

if ( (P%itd_scheme /= 1).and.(P%itd_scheme /= 3) ) then
     call save_sparse(mat_csr, dir_mats // 'temp/' // 'mat_csr_init.dat')
     if (P%itd_scheme == 0) then
     	call dealloc_sparse(mat_csr)
     endif
endif

endif ! endif for the case when baro mat was needed to be assembled

!%======================
!	Internal tide drag	! Done separately bc it requires adding two filled sparse matrices
!%======================
if (P%itd_scheme >= 1) then
	call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
	call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
	call load_alloc_sparse(v2u, dir_mats // 'v2u.dat')
	call load_alloc_sparse(u2v, dir_mats // 'u2v.dat')

	if (P%smooth_type == 0) then
		call load_alloc_vector(H_h, dir_cols // 'H_h.dat')
	else
		call load_alloc_vector(H_h, dir_cols // 'H_sht_h.dat')

		! rough topo for param ITD: raw topo - smooth topo
		if ( (P%itd_scheme == 3).and.(P%baro_on_smoothed/=1) ) then
			call load_alloc_vector(H_sm_sc_h, dir_cols // 'H_h.dat')
		else
			call load_alloc_vector(H_sm_sc_h, dir_cols // 'H_orig_h.dat')
		endif
		H_sm_sc_h = H_sm_sc_h - H_h

	endif

	call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
	call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')

    allocate(dHdph(nu),dHdta(nv), stat = istat)

	call coo_vec_mul(h2uddph, H_h, P%lib, dHdph)
	call coo_vec_mul(h2vddta, H_h, P%lib, dHdta)

	allocate(dHdph_prmtz(nu),dHdta_prmtz(nv), stat = istat)
		if (P%itd_scheme == 3) then
			call coo_vec_mul(h2uddph, H_sm_sc_h, P%lib, dHdph_prmtz)
			call coo_vec_mul(h2vddta, H_sm_sc_h, P%lib, dHdta_prmtz)
		else
			dHdph_prmtz = dHdph
			dHdta_prmtz = dHdta
		endif

	call dealloc_sparse(h2uddph)
	call dealloc_sparse(h2vddta)
	deallocate(H_h)
!	save prior to smoothing
!			call save_vector(dHdph_prmtz, dir_cols // 'dHdph_prmtz.dat')
!			call save_vector(dHdta_prmtz, dir_cols // 'dHdta_prmtz.dat')

! Save dHdph, dHdta for the future use with the IT modelling routine
	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			call save_vector(dHdph*cosh(ta_u)/re, dir_cols // 'gradH_u.dat')
			call save_vector(dHdta*cosh(ta_v)/re, dir_cols // 'gradH_v.dat')
		case (22) ! LAT-LON, TRANSPORT
			call save_vector(dHdph/(re*cos(ta_u)), dir_cols // 'gradH_u.dat')
			call save_vector(dHdta/re, dir_cols // 'gradH_v.dat')
	end select
	deallocate(dHdta,dHdph)

!%======================
!	The residual grad(H) can be very noisy and can introduce noise into the barotropic solution
!	The latter is translated to the IT modelling solver and can be amplified there. Consider smoothing at this step
!%======================
!	IDEA: smooth ||grad(H_residual)|| and |grad_x/grad_y| using running average over itm_avg degrees
!***************************************************
	nf=ceiling(.5*(P%itm_avg*P%nph/360 - 1.)) ! will average over 2*nf+1 cells (extra nf in ph and th directions)
	! generate masks and metrics used below and in itm-solver for moving average
	call metrics_and_masks(nf, P, GD, dir_grid, dir_cols)
	! calculate metrics and masks to later apply running average (u- and v- grids)
	if ( (P%itd_scheme >= 2).and.(P%smooth_itd_param == 1).and.(nf >= 1) ) then

	    call load_alloc_matrix(metrics_ug, dir_grid // 'temp/metrics_ug.dat')
		call load_alloc_matrix(metrics_vg, dir_grid // 'temp/metrics_vg.dat')
		call load_alloc_matrix(mask_ug, dir_grid // 'temp/mask_ug.dat')
		call load_alloc_matrix(mask_vg, dir_grid // 'temp/mask_vg.dat')
		! load up, vp, hp
		call load_alloc_matrix(up, dir_cols // 'up.dat')
		call load_alloc_matrix(vp, dir_cols // 'vp.dat')

		allocate(norm_u(nu), norm_v(nv), ratio_u(nu), ratio_v(nv), stat=istat)
		norm_u=0; norm_v=0; ratio_u=0; ratio_v=0;

		call coo_vec_mul(v2u, dHdta_prmtz, P%lib, norm_u)
		call coo_vec_mul(u2v, dHdph_prmtz, P%lib, norm_v)
		norm_u = (dHdph_prmtz**2 + norm_u**2)**.5
		norm_v = (dHdta_prmtz**2 + norm_v**2)**.5

		where (norm_u>0) ratio_u = abs(dHdph_prmtz)/norm_u
		where (norm_v>0) ratio_v = abs(dHdta_prmtz)/norm_v

		allocate (norm_ug(nph, nth),norm_vg(nph, nth+1), stat=istat)
		norm_ug = 0; norm_vg = 0;
		do cu = 1, nu
		    norm_ug(up(cu,1),up(cu,2)) = norm_u(cu)
		enddo
		do cv = 1, nv
		    norm_vg(vp(cv,1),vp(cv,2)) = norm_v(cv)
		enddo
		allocate (ratio_ug(nph, nth),ratio_vg(nph, nth+1), stat=istat)
		ratio_ug = 0; ratio_vg = 0;
		do cu = 1, nu
		    ratio_ug(up(cu,1),up(cu,2)) = ratio_u(cu)
		enddo
		do cv = 1, nv
		    ratio_vg(vp(cv,1),vp(cv,2)) = ratio_v(cv)
		enddo

		deallocate(norm_u, norm_v, ratio_u, ratio_v)

		! Apply running average to smooth IT drag
!		if (P%sponge_forcing == 0) then
!			call mv_average(Dug, mask_ug*mask_ug_sponge, metrics_ug, nf, P%coor, P%lib)
!			call mv_average(Dvg, mask_vg*mask_vg_sponge, metrics_vg, nf, P%coor, P%lib)
!		else
			call mv_average(norm_ug, mask_ug, metrics_ug, nf, P%coor, P%lib)
			call mv_average(norm_vg, mask_vg, metrics_vg, nf, P%coor, P%lib)

			call mv_average(ratio_ug, mask_ug, metrics_ug, nf, P%coor, P%lib)
			call mv_average(ratio_vg, mask_vg, metrics_vg, nf, P%coor, P%lib)
!		endif
		do cu = 1, nu
		    dHdph_prmtz(cu) = sign( norm_ug(up(cu,1),up(cu,2))*ratio_ug(up(cu,1),up(cu,2)), dHdph_prmtz(cu) )
		enddo
		do cv = 1, nv
		    dHdta_prmtz(cv) = sign( norm_vg(vp(cv,1),vp(cv,2))*ratio_vg(vp(cv,1),vp(cv,2)), dHdta_prmtz(cv) )
		enddo

		deallocate(metrics_ug, metrics_vg, mask_ug, mask_vg, up, vp)
		deallocate(norm_ug, norm_vg, ratio_ug, ratio_vg)

	endif
!	save after smoothing
!			call save_vector(dHdph_prmtz, dir_cols // 'dHdph_prmtz_sm.dat')
!			call save_vector(dHdta_prmtz, dir_cols // 'dHdta_prmtz_sm.dat')

!	allocate(gradH_u_prmtz(nu),gradH_v_prmtz(nv), stat = istat)
	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			dHdph_prmtz = dHdph_prmtz*cosh(ta_u)/re
			dHdta_prmtz = dHdta_prmtz*cosh(ta_v)/re
		case (22) ! LAT-LON, TRANSPORT
			dHdph_prmtz = dHdph_prmtz/(re*cos(ta_u))
			dHdta_prmtz = dHdta_prmtz/re
	end select

	!***********************************************************************************

	call calc_N_on_grid(N_u, N_v, P, nu, nv, nh, N_data_dir, dir_cols, dir_grid, dir_mats)
!		call save_vector( N_u, dir_cols // 'N_u.dat')
!		call save_vector( N_v, dir_cols // 'N_v.dat')


	allocate(Qu(nu), Qv(nv), stat = istat)
	!***********************************************************************************
	select case (cooruv)
	!***********************************************************************************
		case (12) ! MERCATOR, TRANSPORT
			Qu = P%itd_coeff*N_u**2 * (dHdph_prmtz)**2 * cosh(ta_u)
			call init_sparse(DUU, iu, iu, Qu, nu, nu, nu)

		    call right_mult_diag(v2u, (dHdta_prmtz)*cosh(ta_v), nv)
			Qu = P%itd_coeff*N_u**2 * (dHdph_prmtz)
		    call left_mult_diag(v2u, Qu, nu)
			call init_sparse(DVU, v2u%indi, v2u%indj, v2u%vals, v2u%ni, v2u%nj, v2u%nz)

		    call right_mult_diag(u2v, (dHdph_prmtz)*cosh(ta_u), nu)
			Qv = P%itd_coeff*N_v**2 * (dHdta_prmtz)
		    call left_mult_diag(u2v, Qv, nv)
			call init_sparse(DUV, u2v%indi, u2v%indj, u2v%vals, u2v%ni, u2v%nj, u2v%nz)

			Qv = P%itd_coeff*N_v**2 * (dHdta_prmtz)**2 * cosh(ta_v)
			call init_sparse(DVV, iv, iv, Qv, nv, nv, nv)

		case (22) ! LAT-LON, TRANSPORT
			Qu = P%itd_coeff*N_u**2 * (dHdph_prmtz)**2
			call init_sparse(DUU, iu, iu, Qu, nu, nu, nu)

		    call right_mult_diag(v2u, dHdta_prmtz, nv)
			Qu = P%itd_coeff*N_u**2 * (dHdph_prmtz)
		    call left_mult_diag(v2u, Qu, nu)
			call init_sparse(DVU, v2u%indi, v2u%indj, v2u%vals, v2u%ni, v2u%nj, v2u%nz)

		    call right_mult_diag(u2v, dHdph_prmtz, nu)
			Qv = P%itd_coeff*N_v**2 * (dHdta_prmtz)
		    call left_mult_diag(u2v, Qv, nv)
			call init_sparse(DUV, u2v%indi, u2v%indj, u2v%vals, u2v%ni, u2v%nj, u2v%nz)

			Qv = P%itd_coeff*N_v**2 * (dHdta_prmtz)**2
			call init_sparse(DVV, iv, iv, Qv, nv, nv, nv)

	end select

	write(*, '("DUU-DVV, ")', advance = 'no')

!	!***********************************************************************************

	deallocate(Qu, Qv, ta_u, ta_v)
	call deallocate_sparse(v2u)
	call deallocate_sparse(u2v)

	call save_sparse(DUU, dir_mats // 'DUU.dat')
	call save_sparse(DUV, dir_mats // 'DUV.dat')
	call save_sparse(DVU, dir_mats // 'DVU.dat')
	call save_sparse(DVV, dir_mats // 'DVV.dat')

	call deallocate_sparse(DUU)
	call deallocate_sparse(DUV)
	call deallocate_sparse(DVU)
	call deallocate_sparse(DVV)

!==========================================================================================
	if (P%itd_scheme /= 20) then
	!    %=========================================
	!    % again martices for baro tides and not just ITs
	!    %=========================================
		call load_alloc_vector(KU, dir_cols // 'KU.dat')
		call load_alloc_vector(KV, dir_cols // 'KV.dat')

	ncpts=len(cpts)/2

	do ccpt = 1, ncpts

		cpt=cpts(2*ccpt-1:2*ccpt)
		pars = get_pars(cpt)

		call calc_itd_prefactor(ufactor, vfactor, N_u, N_v, dHdph_prmtz, dHdta_prmtz, cpt, P, nu, nv, dir_cols, dir_grid)

		!***********************************************
		!	Add terms into the matrix
		!***********************************************
		call load_alloc_sparse(DUU, dir_mats // 'DUU.dat')
		call load_alloc_sparse(DUV, dir_mats // 'DUV.dat')
		call load_alloc_sparse(DVU, dir_mats // 'DVU.dat')
		call load_alloc_sparse(DVV, dir_mats // 'DVV.dat')

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

		deallocate(ufactor, vfactor)
		call deallocate_sparse(DUU)
		call deallocate_sparse(DUV)
		call deallocate_sparse(DVU)
		call deallocate_sparse(DVV)

		! Convert the tempmat from COO to CSR:
	    call coo2csr(tempmat, tmp_csr, P%lib)
	    call csr_csr_add(mat_csr, real(1., kind=wp), tmp_csr, mkl, mat_csr_cpt)
	    call dealloc_sparse(tmp_csr)
	    call dealloc_sparse(tempmat)
!     Save the COO version of the matrix
!		call csr2coo(mat_csr_cpt, mat, P%lib)
!		call save_sparse(mat, dir_mats // 'temp/' // cpt // '_mat_init.dat')
!		call dealloc_sparse(mat)

!%================================================
!% Store the initial matrix for each cpt in a file
!%================================================
     call save_sparse(mat_csr_cpt, dir_mats // 'temp/' // cpt // '_mat_csr_init.dat')
     call dealloc_sparse(mat_csr_cpt)

	end do

    call dealloc_sparse(mat_csr)
	deallocate( dHdph_prmtz, dHdta_prmtz, N_u, N_v, KU, KV)


	write(*, '("Param ITD (1 mat per cpt), ")', advance = 'no')
!==========================================================================================
	endif ! endif for the case when baro mat was needed to be assembled

end if
!==========================================================================================
!==========================================================================================


!%==================================================
!% define the b.c. diagonal matrix:
!%==================================================
! load up, vp, hp
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')

  allocate(bcdiag(np), stat = istat)
  bcdiag(:) = (/ 1-up(:,3), 1-vp(:,3), (1, j=1,nh) /)

write(*, '("b.c. ")', advance = 'no')

      call save_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
      deallocate(up, vp, bcdiag)
!==========================================================================================

write(*, '("DONE")')

end subroutine
!**********************************************************************************************************
!**********************************************************************************************************
!==========================================================================================
subroutine write_itm_mats(cpts, P, GD, dir_cols, dir_grid, dir_mats)

!     character(len = *) :: N_data_dir
     type(params), intent(in)	:: P
     type(grid_dims), intent(in):: GD
     character(len=*) :: dir_cols, dir_grid, dir_mats

     integer :: nu, nv, nh, np, cooruv!, nph, nth
     real(wp)	:: re

     integer, target  :: iu(GD%nu), iv(GD%nv), ih(GD%nh)
     real(wp), allocatable ::  ta_u(:), ta_v(:), ta_h(:)!, ph_vg(:), ta_ug(:)

!     real    ::      T1, T2
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     integer ::      j, istat

     character(len=*) :: cpts
!     character(len=2) :: cpt
!	 type (tide_params) :: pars
     integer :: nmodes, cmode !, ncpts, ccpt


!     real(wp), allocatable :: H_u(:), H_v(:), H_h(:), KU(:), KV(:)
!     real(wp), allocatable :: dHdph(:), dHdta(:), N_u(:), N_v(:)

!     real(wp), allocatable	:: issponge_u(:, :), issponge_v(:, :), issponge_h(:, :)
	 real(dp), allocatable	:: c2n_v(:,:), c2n_u(:,:)

	 type (triplet)			:: tempmat
     integer				:: nmat, itm_filled, itm_filled_0 ! points how many elements out of nmat are already filled
!	 integer, allocatable	:: ufactor(:), vfactor(:)

     type (triplet) :: u2vf, v2uf, h2uddph, h2vddta, u2hddph, v2hddta

     type (triplet) :: itm_mat
	 type (csr) :: itm_mat_csr

     type (tide_params)	   	:: t_params
     integer				:: tide_types(P%ncpts), ntypes, ctypes

!%=============================================
! shortcuts
!%=============================================
!g = P%g
re = P%re
nmodes = P%n_modes
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
!%=============================================

!    %==========================
!    % set up the basic matrices
!    %==========================
!	Should already be done when called baro_uvhT_mat(nu, nv, nh, P, nmat, dir_cols, dir_mats)


     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (j), j=1,nv ) /)
     ih = (/ ( (j), j=1,nh ) /)

    call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
    call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')
	call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')

!    %=====================================
!    % rotation, c.o.mass and RHS (forcing)
!    %=====================================
!%=============================================
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
! Create the linear IT model matrix:
!	Target: 1) rotation		 2) pressure gradients		3) divergence terms		  ( 4) diagonal i*omega & sponge damping terms ) -- freq dependent, added in the mail loop
!	nmat = v2uf%nz + u2vf%nz + h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz + nu + nv

	call load_alloc_sparse(v2uf, dir_mats // 'v2uf.dat')
	call load_alloc_sparse(u2vf, dir_mats // 'u2vf.dat')
	call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
	call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
	call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
	call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')

!	nmat = v2uf%nz + u2vf%nz + h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz + nu + nv + nh
	nmat = v2uf%nz + u2vf%nz + h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz ! ( + nu + nv + nh -- diag terms are added later)
	np = nu + nv + nh
	itm_filled_0 = 0

	call dealloc_sparse(u2hddph)
	call dealloc_sparse(v2hddta)

	call init_sparse_0(itm_mat, np,np, nmat)

write(*, '("Assemble global IT gen matrix: 1) common part for all modes, ")', advance = 'no')
!%===============
!	Rotation
!%===============
!	1) v2uf
	v2uf%vals = - v2uf%vals ! multiply by -1
	call concatenate_sparse(itm_mat, v2uf, 0, nu, itm_filled_0)
!	2) u2vf
    call concatenate_sparse(itm_mat, u2vf, nu, 0, itm_filled_0)

    call dealloc_sparse(v2uf)
	call dealloc_sparse(u2vf)
!%==============================
!	Pressure gradients terms
!%==============================
		select case (cooruv)
			case (12) ! MERCATOR, TRANSPORT
			!	1) Grad P on u-grid
				h2uddph%vals = 1/re*h2uddph%vals
			!	2) Grad P on v-grid
				h2vddta%vals = 1/re*h2vddta%vals
			case (22) ! LAT-LON, TRANSPORT
			!	1) Grad P on u-grid
		    	call left_mult_diag(h2uddph, 1/re/cos(ta_u), nu)
			!	2) Grad P on v-grid
				h2vddta%vals = 1/re*h2vddta%vals
		end select
		call concatenate_sparse(itm_mat, h2uddph, 0, nu+nv, itm_filled_0)
		call concatenate_sparse(itm_mat, h2vddta, nu, nu+nv, itm_filled_0)

		call dealloc_sparse(h2uddph)
		call dealloc_sparse(h2vddta)

call save_sparse(itm_mat, dir_mats // 'itm/temp/itm_mat_0.dat')
call dealloc_sparse(itm_mat)

!%===================================
! Advection  (+ linear dumping at the sponge) -- in the main cycle
!%===================================
write(*, '("2) individual (1 mat per mode and cpt_2)... ")', advance = 'no')

!call load_alloc_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//cpt(2:2)//'.dat')
!call load_alloc_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//cpt(2:2)//'.dat')
!call load_alloc_matrix(issponge_h, dir_grid//'itm/'//'issponge_h'//cpt(2:2)//'.dat')

!	Generated TWO different cn and sponges for semidiurnal and diurnal constituents
	! check which are present
	tide_types = 0
	do j = 1, P%ncpts
		if ( cpts(2*j:2*j) .eq. '1') then ! diurnal
			tide_types(j) = 1
		elseif ( cpts(2*j:2*j) .eq. '2') then ! semidiurnal
			tide_types(j) = 2
		endif
	enddo
	call unique(tide_types, ntypes)

do ctypes = 1, ntypes ! a tiny cycle to make separate modified cn and sponge for diurnal and semidiurnal tides
	call load_alloc_matrix(c2n_u, dir_grid //'itm/' // 'c2n_u'//tostring(tide_types(ctypes))//'.dat')
	call load_alloc_matrix(c2n_v, dir_grid //'itm/' // 'c2n_v'//tostring(tide_types(ctypes))//'.dat')

	do cmode = 1, nmodes
		write(*, '(i2,a) ', advance='no') cmode, ', '

		call load_alloc_sparse(itm_mat, dir_mats // 'itm/temp/itm_mat_0.dat')
		itm_filled = itm_filled_0

		!%===================================
		!	Advection gradients (c.o. mass)
		!%===================================
			call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
			call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')
			select case (cooruv)
				case (12) ! MERCATOR, TRANSPORT
					!	1) d/dx( c2n * u )
			    		call left_mult_diag(u2hddph, (1/re)*cosh(ta_h)**2, nh)
			    		call right_mult_diag(u2hddph, c2n_u(:,cmode), nu)
					!	2) d/dy( c2n * v )
						call left_mult_diag(v2hddta, (1/re)*cosh(ta_h)**2, nh)
						call right_mult_diag(v2hddta, c2n_v(:,cmode), nv)
				case (22) ! LAT-LON, TRANSPORT
					!	1) d/dx( c2n * u )
			    		call left_mult_diag(u2hddph, (1/re)*1/cos(ta_h), nh)
			    		call right_mult_diag(u2hddph, c2n_u(:,cmode), nu)
					!	2) d/dy( c2n * v )
						call left_mult_diag(v2hddta, (1/re)*1/cos(ta_h), nh)
						call right_mult_diag(v2hddta, c2n_v(:,cmode)*cos(ta_v), nv)
			end select
			call concatenate_sparse(itm_mat, u2hddph, nu+nv, 0, itm_filled)
			call concatenate_sparse(itm_mat, v2hddta, nu+nv, nu, itm_filled)

			call dealloc_sparse(u2hddph)
			call dealloc_sparse(v2hddta)

		!%===================================
		!	Diagonal sponge damping terms ( -- in the main cycle)
		!%===================================
		!--------------- IMPORTANT ---------------
		! i*omega & sponge ARE included together at a later stage as a diag mat (-i*omega + omega*P%sponge_damp)

	!		call init_sparse(tempmat, iu, iu, P%sponge_damp*issponge_u(:, cmode), nu, nu, nu)
	!		call concatenate_sparse(itm_mat, tempmat, 0, 0, itm_filled)
	!		call dealloc_sparse(tempmat)
	!
	!		call init_sparse(tempmat, iv, iv, P%sponge_damp*issponge_v(:, cmode), nv,nv, nv)
	!		call concatenate_sparse(itm_mat, tempmat, nu, nu, itm_filled)
	!		call dealloc_sparse(tempmat)

	!!		write(*, '("sponge, ")', advance = 'no')

	!%================================================
	!% Store the initial matrix for each cpt in a file
	!%================================================
	! Convert mat from COO to CSR:
	     call coo2csr(itm_mat, itm_mat_csr, P%lib)
	!     call save_sparse(mat, dir_mats // 'temp/mat_init.dat')
	     call dealloc_sparse(itm_mat)

	     call save_sparse(itm_mat_csr, dir_mats // 'itm/temp/' // 'itm_mat_csr_init_'// &
	     										  tostring(tide_types(ctypes))//'_'//tostring(cmode)//'.dat')
	     call dealloc_sparse(itm_mat_csr)

	enddo
	deallocate(c2n_u, c2n_v)
enddo
!==========================================================================================
write(*, '("DONE")')


end subroutine

!**********************************************************************************************************
!**********************************************************************************************************

subroutine baro_uvhT_mat(nu, nv, nh, P, nmat, dir_cols, dir_mats)

    implicit none

!% Produces a linear inversion/propagator matrix for a
!% column uvh using the TRANSPORT VELOCITY formulation.
!% Includes LINEAR friction
!% DOES NOT include scalar SAL
!% No frequency dependent d/dt or internal tide drag terms.
     character(len=*) :: dir_cols, dir_mats
!     character(len = *) :: N_data_dir
     type(params) :: P
     integer, intent(in) :: nu, nv, nh
     integer, intent(out) :: nmat

	integer		:: cooruv !differentiate between lat-lon/MERC, velocity/transport formulations
!     integer, allocatable :: up(:, :), vp(:, :)
     real(wp)			  :: g, re
     real(wp), allocatable ::  ta_u(:), ta_v(:), ta_h(:), ones(:)
     real(wp), allocatable :: H_u(:), H_v(:)!, H_h(:)
     integer	:: istat

     type (triplet) :: h2uddph, h2vddta, u2hddph, v2hddta
!     type (triplet) :: GHU, GHV, JUH, JVH, DUU, DVU, DUV, DVV
!     real(wp), allocatable :: KU(:), KV(:)

!%=============================================
! shortcuts
!%=============================================
g = P%g
re = P%re
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
!%=============================================



!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
!	nmat = h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
!	realized in steps 1-2 and 5-6
	nmat = nu + nv

	call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
	call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')
	call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')

select case (cooruv)
!***********************************************************************************
	case (12) ! MERCATOR, TRANSPORT
!***********************************************************************************
	!%================================================
	!	Pressure gradients (gravity like potentials)
	!%================================================
	!	1) GHU
		call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
		call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	    call left_mult_diag(h2uddph, (g/re)*H_u, nu)

		call save_sparse(h2uddph, dir_mats // 'GHU.dat')
		nmat = nmat + h2uddph%nz
		call dealloc_sparse(h2uddph)
		deallocate(H_u)

		write(*, '("GHU, ")', advance = 'no')

	!	2) GHV
		call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
		call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
	    call left_mult_diag(h2vddta, (g/re)*H_v, nv)

		call save_sparse(h2vddta, dir_mats // 'GHV.dat')
		nmat = nmat + h2vddta%nz
		call dealloc_sparse(h2vddta)
		deallocate(H_v)

		write(*, '("GHV, ")', advance = 'no')

	!%================================================
	!	Metrics multipliers
	!%================================================
	!	3) KU
		call save_vector(1/cosh(ta_u), dir_cols // 'KU.dat')
		write(*, '("KU, ")', advance = 'no')
	!	4) KV
		call save_vector(1/cosh(ta_v), dir_cols // 'KV.dat')
		write(*, '("KV, ")', advance = 'no')

	!%================================================
	!	Advection gradients (c.o. mass)
	!%================================================
	!	5) JUH
		call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
	    call left_mult_diag(u2hddph, (1/re)*cosh(ta_h)**2, nh)

		call save_sparse(u2hddph, dir_mats // 'JUH.dat')
		nmat = nmat + u2hddph%nz
		call dealloc_sparse(u2hddph)

	!	6) JVH
		call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')
		call left_mult_diag(v2hddta, (1/re)*cosh(ta_h)**2, nh)

		call save_sparse(v2hddta, dir_mats // 'JVH.dat')
		nmat = nmat + v2hddta%nz
		call dealloc_sparse(v2hddta)

!***********************************************************************************
	case (22) ! LAT-LON, TRANSPORT
!***********************************************************************************
	!%================================================
	!	Pressure gradients (gravity like potentials)
	!%================================================
	!	1) GHU
		call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
		call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	    call left_mult_diag(h2uddph, (g/re)*H_u/cos(ta_u), nu)

		call save_sparse(h2uddph, dir_mats // 'GHU.dat')
		nmat = nmat + h2uddph%nz
		call dealloc_sparse(h2uddph)
		deallocate(H_u)

		write(*, '("GHU, ")', advance = 'no')

	!	2) GHV
		call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
		call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
	    call left_mult_diag(h2vddta, (g/re)*H_v, nv)

		call save_sparse(h2vddta, dir_mats // 'GHV.dat')
		nmat = nmat + h2vddta%nz
		call dealloc_sparse(h2vddta)
		deallocate(H_v)

		write(*, '("GHV, ")', advance = 'no')

	!%================================================
	!	Metrics multipliers
	!%================================================
	!	3) KU
		allocate(ones(nu), stat = istat)
		ones = 1
		call save_vector(ones, dir_cols // 'KU.dat')
		deallocate(ones)
		write(*, '("KU, ")', advance = 'no')
	!	4) KV
		allocate(ones(nv), stat = istat)
		ones = 1
		call save_vector(ones, dir_cols // 'KV.dat')
		deallocate(ones)
		write(*, '("KV, ")', advance = 'no')

	!%================================================
	!	Advection gradients (c.o. mass)
	!%================================================
	!	5) JUH
		call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
	    call left_mult_diag(u2hddph, (1/re)*1/cos(ta_h), nh)

		call save_sparse(u2hddph, dir_mats // 'JUH.dat')
		nmat = nmat + u2hddph%nz
		call dealloc_sparse(u2hddph)

	!	6) JVH
		call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')
		call left_mult_diag(v2hddta, (1/re)*1/cos(ta_h), nh)
		call right_mult_diag(v2hddta, cos(ta_v), nv)

		call save_sparse(v2hddta, dir_mats // 'JVH.dat')
		nmat = nmat + v2hddta%nz
		call dealloc_sparse(v2hddta)

end select

	deallocate(ta_h, ta_u, ta_v)

end subroutine baro_uvhT_mat

!==========================================================================================
!==========================================================================================

subroutine generate_global_matrices(dir_grid, dir_cols, dir_mats, P, GD)

     character(len=*) :: dir_grid, dir_cols, dir_mats
     type(params) :: P
     type(grid_dims):: GD

     logical :: flag_itd_smooth, flag_baro_smooth

!     real(dp), allocatable :: H_hg(:, :), H_sht_hg(:, :)


!%=======================================
!% write H, ta, heq etc. on the columns:
!%=======================================
	 flag_itd_smooth = (P%itd_scheme>0).and.(P%smooth_type>0).and.(P%sht_smooth_H >= 0)
	 flag_baro_smooth = (P%smooth_type>0).and.(P%baro_on_smoothed == 1)
     call write_global_ccols_col(GD%nph, GD%nu, GD%nv, GD%nh, dir_grid, dir_cols, flag_itd_smooth, flag_baro_smooth)

!     deallocate(H_hg,H_sht_hg) ! don't need it anymore, free space
!%====================================
!% generate and write sparse matrices
!%====================================
     call write_global_cmats_col(GD%nph, GD%nu, GD%nv, GD%nh, P%latP, P%omega, dir_grid, dir_cols, dir_mats)

!     deallocate(ta_ug, ta_vg, ph_ug, ph_vg, H_hg)
!     deallocate(up, vp, hp, iu_ug, iv_vg, ih_hg)

     write(*, '("done.")')

end subroutine

!    %========================================================
!    % BACIS matrices
!    %========================================================
!==========================================================================================
subroutine write_global_cmats_col(nph, nu, nv, nh, latP, omega, dir_grid, dir_cols, dir_mats)
!%===============================================================
!% load the grid: H_hg, ta_ug, ta_vg, ph_ug, ph_vg, th_ug, th_vg,
!%                iu_ug, iv_vg, ih_hg
!%===============================================================
    implicit none

!% Generate some sparse matrices for operations on an Arakawa C grid.
!%
!%   h2uddph, h2u
!%   h2vddta, h2v
!%   u2hddph, u2h
!%   v2hddta, v2h
!%   u2v, u2vf, u2vfsp
!%   v2u, v2uf, v2ufsp

     integer                   ::     nph!, nth
     real(dp), allocatable     ::     th_vg(:), ta_vg(:), ph_ug(:)

     integer, allocatable :: up(:, :), vp(:, :), hp(:, :)
     integer, intent(in)  ::     nu, nv, nh!, np
     integer              :: nmax
     integer, allocatable :: iu_ug(:, :), iv_vg(:, :), ih_hg(:, :)
     real(dp), allocatable :: H_u(:), H_v(:)!, H_h(:)


     integer, pointer :: ivals(:), jvals(:), ispvals(:), jspvals(:)
     real(dp), pointer :: vals(:), vals1(:), fvals(:), fspvals(:)

     real(dp)  ::      dph, dta
     integer ::      istat, statusj, status1, status2
     integer ::      cu, cv, ch, cth, cph
     integer ::      cphl, cphr, ctha, cthb, nvs
     integer ::      hpl, hpr, hpa, hpb
     integer ::      upl, upr
     integer ::      vpa, vpb
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     real(dp), intent(in) :: omega, latP
     real(dp)			  :: th0

     character(len = *) :: dir_grid, dir_cols, dir_mats

!     real(dp), dimension(4)   :: thc, phc

     type (triplet_dp) :: h2uddph, h2u, h2vddta, h2v, u2hddph, u2h, v2hddta, v2h
     type (triplet_dp) :: u2v, u2vf, u2vfsp, v2u, v2uf, v2ufsp

!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
! load up, vp, hp
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

     call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat') ! in Mercator coords it is ta
     														! in lat-lon coords it is th
     call load_alloc_vector(th_vg, dir_grid // 'th_vg.dat')
! load iu_ug, iv_vg, ih_hg
     call load_alloc_matrix(iu_ug, dir_grid // 'iu_ug.dat')
     call load_alloc_matrix(iv_vg, dir_grid // 'iv_vg.dat')
     call load_alloc_matrix(ih_hg, dir_grid // 'ih_hg.dat')
! Load columns H_u, H_v
     call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
     call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
!%=================
!% make shortcuts:
!%=================

th0 = d2r(latP)
!ph0 = d2r(lonP)

!dph=2*pi/nph;    ! the distance between ph gridpoints
dph = ph_ug(2)-ph_ug(1)
!dta=dph;    ! the distance between ta gridpoints
dta = ta_vg(2)-ta_vg(1)

!==========================================================================================
! DISP MODULE SETTINGS
!==========================================================================================
call tostring_set(rfmt='F12.1')

!     SECOND ORDER FINITE DIFFERENCE SCHEME

write(*, '("Making matrices:")')

nmax = max(nu, nv, nh)
    if (associated(ivals)) deallocate(ivals) ! 4*nmax for 4th order scheme
    allocate(ivals(2*nmax), stat = istat)
    if (associated(jvals)) deallocate(jvals)
    allocate(jvals(2*nmax), stat = statusj)
    if (associated(vals)) deallocate(vals)
    allocate(vals(2*nmax), stat = status1)
    if (associated(vals1)) deallocate(vals1)
    allocate(vals1(2*nmax), stat = status2)

    ivals = 0
    jvals = 0
    vals = 0
    vals1 = 0

    nvs = 0


write(*, '(" - d/dph and 1 for h-grid functions, evaluating on the u-grid ")', advance = 'no')

!%=============
!% initialize:
!%=============
!!!!!!!!!!!!!!       h2uddph, h2u
!%=============

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cu = 1,nu

       cph=up(cu,1)
       cth=up(cu,2)

      cphr=cph
      cphl=1+modulo((cph-1)-1,nph)
       hpl=ih_hg(cphl,cth)
       hpr=ih_hg(cphr,cth)

       if ((hpl /= 0).and.(hpr /= 0)) then !% we are not at a Western/Eastern boundary

      !% 4th order FD
!              cphll=1+mod((cphl-1)-1,nph)
!              cphrr=1+mod((cphr+1)-1,nph);
!
!              hpll=ih_hg(cth,cphll); hprr=ih_hg(cth,cphrr);
!              if hpll > 0 & hprr > 0 & fdflag == 4
!                nvs=nvs+4;
!                ivals(nvs-3:nvs)=cu;
!                jvals(nvs-3:nvs)=[hpll hpl hpr hprr];
!                vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dph;
!                vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!              else
                nvs=nvs+2
                ivals(nvs-1:nvs)=cu
                jvals(nvs-1:nvs)=(/hpl, hpr/)
                vals(nvs-1:nvs)=(/-1., 1./)/dph
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
       endif

    enddo

call init_sparse_dp(h2uddph,ivals,jvals,vals,nu,nh,nvs)
call init_sparse_dp(h2u,ivals,jvals,vals1,nu,nh,nvs)

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       h2vddta, h2v
!%=============

write(*, '(" - d/dta and 1 for h grid functions, evaluating on the v-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cv = 1,nv

       cph=vp(cv,1)
       cth=vp(cv,2)

      ctha=cth
      cthb=cth-1
       hpa=ih_hg(cph,ctha)
       hpb=ih_hg(cph,cthb)

       if ((hpa /= 0).and.(hpb /= 0)) then !% we are not at a Northern/Southern boundary

      !% 4th order FD
!         cthaa=ctha+1; cthbb=cthb-1;
!         hpaa=ih_hg(cthaa,cph); hpbb=ih_hg(cthbb,cph);
!         if hpaa > 0 && hpbb > 0 && fdflag == 4
!           nvs=nvs+4;
!           ivals(nvs-3:nvs)=cv;
!           jvals(nvs-3:nvs)=[hpbb hpb hpa hpaa];
!           vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dta;
!           vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!         else
                nvs=nvs+2
                ivals(nvs-1:nvs)=cv
                jvals(nvs-1:nvs)=(/hpb, hpa/)
                vals(nvs-1:nvs)=(/-1., 1./)/dta
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
       endif

    enddo

call init_sparse_dp(h2vddta,ivals,jvals,vals,nv,nh,nvs)
call init_sparse_dp(h2v,ivals,jvals,vals1,nv,nh,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')


!%=============
!!!!!!!!!!!!!!       u2hddph, u2h
!%=============

write(*, '(" - d/dph and 1 for u-grid functions, evaluating on the h-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do ch = 1,nh

       cph=hp(ch,1)
       cth=hp(ch,2)

      cphl=cph
      cphr=1+modulo((cph+1)-1,nph)
       upl=iu_ug(cphl,cth)
       upr=iu_ug(cphr,cth)


!       if ((hpa /= 0).and.(hpb /= 0)) then !% we are not at a Northern/Southern boundary

      !% 4th order FD
!      cphll=1+mod((cphl-1)-1,nph); cphrr=1+mod((cphr+1)-1,nph);
!      upll=iu_ug(cth,cphll); uprr=iu_ug(cth,cphrr);
!       if upll > 0 & uprr > 0 & fdflag == 4
!         nvs=nvs+4;
!         ivals(nvs-3:nvs)=ch;
!         jvals(nvs-3:nvs)=[upll upl upr uprr];
!         vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dph;
!         vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!       else
                nvs=nvs+2
                ivals(nvs-1:nvs)=ch
                jvals(nvs-1:nvs)=(/upl, upr/)
                vals(nvs-1:nvs)=(/-1., 1./)/dph
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
!       endif

    enddo

call init_sparse_dp(u2hddph,ivals,jvals,vals,nh,nu,nvs)
call init_sparse_dp(u2h,ivals,jvals,vals1,nh,nu,nvs)



call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       v2hddta, v2h
!%=============

write(*, '(" - d/dta and 1 for v-grid functions, evaluating on the h-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do ch = 1,nh

       cph=hp(ch,1)
       cth=hp(ch,2)

      ctha=cth+1
      cthb=cth
       vpa=iv_vg(cph,ctha)
       vpb=iv_vg(cph,cthb)

!       cthaa=ctha+1; cthbb=cthb-1;
!       vpaa=iv_vg(cthaa,cph); vpbb=iv_vg(cthbb,cph);
!       if vpaa > 0 & vpbb > 0 & fdflag == 4
!         nvs=nvs+4;
!         ivals(nvs-3:nvs)=ch;
!         jvals(nvs-3:nvs)=[vpbb vpb vpa vpaa];
!         vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dta;
!         vals1(nvs-3:nvs)=[-1 9 9 -1]/18;
!       else
                nvs=nvs+2
                ivals(nvs-1:nvs)=ch
                jvals(nvs-1:nvs)=(/vpb, vpa/)
                vals(nvs-1:nvs)=(/-1., 1./)/dta
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
!       endif

    enddo

call init_sparse_dp(v2hddta,ivals,jvals,vals,nh,nv,nvs)
call init_sparse_dp(v2h,ivals,jvals,vals1,nh,nv,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!**********************************************************************************************************
!**********************************************************************************************************
    deallocate(ivals, jvals, vals, vals1)
    allocate(ivals(4*nmax), jvals(4*nmax), vals(4*nmax), stat = status1)

    if (associated(fvals)) deallocate(fvals)
    allocate(fvals(4*nmax), stat = status2)

    if (associated(ispvals)) deallocate(ispvals) ! nmax for....
    allocate(ispvals(2*nmax), stat = istat)
    if (associated(jspvals)) deallocate(jspvals)
    allocate(jspvals(2*nmax), stat = statusj)
    if (associated(fspvals)) deallocate(fspvals)
    allocate(fspvals(2*nmax), stat = status1)

    ivals = 0
    jvals = 0
    vals = 0

    ispvals = 0
    jspvals = 0
    fspvals = 0
!**********************************************************************************************************
!**********************************************************************************************************
!%=============
!!!!!!!!!!!!!!       v2u, v2uf, v2ufsp
!%=============

write(*, '(" - v & fv,     			     evaluating on the u-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cu = 1,nu

       cph=up(cu,1)
       cth=up(cu,2)

       if (up(cu,3) == 0) then

       ctha=cth+1
       cthb=cth
       cphr=cph
       cphl=1+modulo((cph-1)-1,nph)

                nvs=nvs+4
                ivals(nvs-3:nvs)=cu
                jvals(nvs-3:nvs)=(/iv_vg(cphr,ctha), iv_vg(cphr,cthb), iv_vg(cphl,ctha), iv_vg(cphl,cthb)/)
                vals(nvs-3:nvs)=(/0.25, 0.25, 0.25, 0.25/)

                fvals(nvs-3:nvs)=0.25*calc_f_rot((/th_vg(ctha), th_vg(cthb), th_vg(ctha), th_vg(cthb)/), &
                                                 (/ph_ug(cph), ph_ug(cph), ph_ug(cph), ph_ug(cph)/), &
                                                 th0, omega)

!		    if flags.baro.num.f(2) == 2
!		      cf=1;
!		      while vp(jvals(nvs-4+cf)) == 1
!		        cf=cf+1;
!		      end
!		      ispvals(nvs/4)=cu;
!		      jspvals(nvs/4)=jvals(nvs-4+cf);
!		      fspvals(nvs/4)=4*fvals(nvs-4+cf);
!    elseif flags.baro.num.f(2) == 3
      ispvals(nvs/2-1:nvs/2)=cu;
      jspvals(nvs/2-1:nvs/2)=(/jvals(nvs-3), jvals(nvs)/)
      fspvals(nvs/2-1:nvs/2)=2*(/fvals(nvs-3), fvals(nvs)/)
!    end
       endif

    enddo

call init_sparse_dp(v2u,ivals,jvals,vals,nu,nv,nvs)
call init_sparse_dp(v2uf,ivals,jvals,fvals,nu,nv,nvs)

!if flags.baro.num.f(2) == 2
!  ispvals=ispvals(1:nvs/4);
!  jspvals=jspvals(1:nvs/4);
!  fspvals=fspvals(1:nvs/4);
!  v2ufsp=sparse(ispvals,jspvals,fspvals,nu,nv);
!  save(matfile,'-append','v2ufsp');
!  clear v2ufsp
!elseif flags.baro.num.f(2) == 3

call init_sparse_dp(v2ufsp,ispvals,jspvals,fspvals,nu,nv,nvs/2) ! nvs/2 elements!

!  save(matfile,'-append','v2ufsp');
!  clear v2ufsp
!end


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       u2v, u2vf, u2vfsp
!%=============

write(*, '(" - u & fu,     			     evaluating on the v-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cv = 1,nv

       cph=vp(cv,1)
       cth=vp(cv,2)

       if (vp(cv,3) == 0) then

       ctha=cth
       cthb=cth-1
       cphr=1+modulo((cph+1)-1,nph)
       cphl=cph


                nvs=nvs+4
                ivals(nvs-3:nvs)=cv
                jvals(nvs-3:nvs)=(/iu_ug(cphr,ctha), iu_ug(cphr,cthb), iu_ug(cphl,ctha), iu_ug(cphl,cthb)/)
                vals(nvs-3:nvs)=(/0.25, 0.25, 0.25, 0.25/)
                fvals(nvs-3:nvs)=0.25*calc_f_rot((/th_vg(cth), th_vg(cth), th_vg(cth), th_vg(cth)/), &
                                                 (/ph_ug(cphr), ph_ug(cphr), ph_ug(cphl), ph_ug(cphl)/), &
                                                 th0, omega)

!    if flags.baro.num.f(2) == 2
!      cf=1;
!      while up(jvals(nvs-4+cf)) == 1
!        cf=cf+1;
!      end
!      ispvals(nvs/4)=cv;
!      jspvals(nvs/4)=jvals(nvs-4+cf);
!      fspvals(nvs/4)=4*fvals(nvs-4+cf);
!              elseif flags.baro.num.f(2) == 3
                ispvals(nvs/2-1:nvs/2)=cv;
                jspvals(nvs/2-1:nvs/2)=(/jvals(nvs-3), jvals(nvs)/)
                fspvals(nvs/2-1:nvs/2)=2*(/fvals(nvs-3), fvals(nvs)/)
!              end
       endif

    enddo

call init_sparse_dp(u2v,ivals,jvals,vals,nv,nu,nvs)
call init_sparse_dp(u2vf,ivals,jvals,fvals,nv,nu,nvs)

!if flags.baro.num.f(2) == 2
!  ispvals=ispvals(1:nvs/4);
!  jspvals=jspvals(1:nvs/4);
!  fspvals=fspvals(1:nvs/4);
!  u2vfsp=sparse(ispvals,jspvals,fspvals,nv,nu);
!  save(matfile,'-append','u2vfsp');
!  clear u2vfsp
!elseif flags.baro.num.f(2) == 3

call init_sparse_dp(u2vfsp,ispvals,jspvals,fspvals,nv,nu,nvs/2) ! nvs/2 elements!

!  save(matfile,'-append','u2vfsp');
!  clear u2vfsp
!end


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%================================================
!% Deallocate what's possible
!%================================================
      deallocate(iu_ug, iv_vg, ih_hg, ph_ug, th_vg, ta_vg)

!%===================================================================================================================
!% save sparse matrices h2uddph, h2u, h2vddta, h2v, u2hddph, u2h, v2hddta, v2h, u2v, u2vf, u2vfsp, v2u, v2uf, v2ufsp
!%===================================================================================================================
write(*, '("Saving matrices: ")', advance = 'no')
      call save_cmat_dp(h2uddph, dir_mats // 'h2uddph.dat')
      call save_cmat_dp(h2u, dir_mats // 'h2u.dat')
      call save_cmat_dp(h2vddta, dir_mats // 'h2vddta.dat')
      call save_cmat_dp(h2v, dir_mats // 'h2v.dat')
      call save_cmat_dp(u2hddph, dir_mats // 'u2hddph.dat')
      call save_cmat_dp(u2h, dir_mats // 'u2h.dat')
      call save_cmat_dp(v2hddta, dir_mats // 'v2hddta.dat')
      call save_cmat_dp(v2h, dir_mats // 'v2h.dat')
      call save_cmat_dp(u2v, dir_mats // 'u2v.dat')
      call save_cmat_dp(v2u, dir_mats // 'v2u.dat')

      call save_cmat_dp(u2vfsp, dir_mats // 'u2vfsp.dat')
      call save_cmat_dp(v2ufsp, dir_mats // 'v2ufsp.dat')

!%================================================
!% Adjust the Coriolis matrices to conserve energy
!%================================================
      call write_modified_fmats_col(u2vf, v2uf, H_u, H_v, nu, nv)

      call save_cmat_dp(u2vf, dir_mats // 'u2vf.dat')
      call save_cmat_dp(v2uf, dir_mats // 'v2uf.dat')

call tostring_set_factory()

end subroutine write_global_cmats_col

!==========================================================================================
!==========================================================================================

subroutine write_modified_fmats_col(u2vf, v2uf, H_u, H_v, nu, nv)
!%====================================================================
!% Different when solving in velocity or volume transport formulations
!%====================================================================
    implicit none

     integer, intent(in)   ::     nu, nv
     real(dp), dimension(:) :: H_u(nu), H_v(nv)

     type (triplet_dp)   :: u2vf, v2uf

!     integer ::      cu, cv, cz

!	Formulation in terms of velocity
!if flags.baro.num.uvform(domflag) == 1
!  v2uf=spdiags(1./sqrt(H_u),0,nu,nu)*v2uf*spdiags(sqrt(H_v),0,nv,nv);
!  u2vf=spdiags(1./sqrt(H_v),0,nv,nv)*u2vf*spdiags(sqrt(H_u),0,nu,nu);
!elseif flags.baro.num.uvform(domflag) == 2
!	Formulation in terms of volume transport

!  v2uf=spdiags(sqrt(H_u),0,nu,nu)*v2uf*spdiags(1./sqrt(H_v),0,nv,nv)
   call left_right_mult_diag_dp(v2uf, (H_u)**0.5, nu, 1/((H_v)**0.5), nv)

!  u2vf=spdiags(sqrt(H_v),0,nv,nv)*u2vf*spdiags(1./sqrt(H_u),0,nu,nu)
   call left_right_mult_diag_dp(u2vf, (H_v)**0.5, nv, 1/((H_u)**0.5), nu)

!end

end subroutine write_modified_fmats_col

!==========================================================================================
!==========================================================================================

function calc_f_rot(thc,phc,th0,omega)

implicit none

real(dp), intent(in)           :: thc(:), phc(:), th0, omega
real(dp), dimension(size(phc)) :: calc_f_rot

!% calculates the Coriolis parameters at a point on the
!% computational lat-lon grid. The computational grid is
!% shifted by ph0 in lon, and dropped down th0 from the
!% pole (in the direction of ph0).

calc_f_rot = 2*omega*(-sin(th0)*cos(thc)*cos(phc)+cos(th0)*sin(thc))

end function calc_f_rot

!**********************************************************************************************************
!**********************************************************************************************************

subroutine write_global_ccols_col(nph, nu, nv, nh, dir_grid, dir_cols, flag_itd_smooth, flag_baro_smooth)

    implicit none

!% Generate some useful vectors when using an Arakawa C grid:
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   H_sht_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points

     integer, intent(in) ::     nph!, nth
     integer, intent(in) ::     nu, nv, nh!, np
     real(dp), allocatable    ::     H_hg(:, :), H_sht_hg(:, :)
     real(dp), allocatable  :: ta_ug(:), ta_vg(:)!, ph_vg(:)
     integer, allocatable   :: up(:, :), vp(:, :), hp(:, :)

     real(dp), pointer :: ta_u(:), ta_v(:), ta_h(:)
     real(dp), pointer :: H_u(:), H_v(:), H_h(:), H_sht_u(:), H_sht_v(:), H_sht_h(:)

     integer ::      statusu, statusv, statush
     integer ::      cu, cv, ch, cth, cph, cphl, cphr, ctha, cthb
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     character(len=*) :: dir_grid, dir_cols
     logical 		:: flag_itd_smooth, flag_baro_smooth
!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

	 call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
	!	SHT-smoothed bottom topography
	if ( flag_itd_smooth ) then
	    call load_alloc_matrix(H_sht_hg, dir_grid // 'H_sht_hg.dat')

	    if (associated(H_sht_u)) deallocate(H_sht_u)
    	allocate(H_sht_u(nu), stat = statusu)
		if (associated(H_sht_v)) deallocate(H_sht_v)
    	allocate(H_sht_v(nv), stat = statusv)
    	if (associated(H_sht_h)) deallocate(H_sht_h)
    	allocate(H_sht_h(nh), stat = statush)
	endif


!write(*, '("Making vectors: ")')
write(*, '("Making vectors: ")', advance = 'no')

!write(*, '(" - tau,          evaluating on the grid-points...... ")', advance = 'no')
write(*, '("tau, ")', advance = 'no')
!%=============
!% initialize:
!%=============

!     tau
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(ta_u)) deallocate(ta_u)
    allocate(ta_u(nu), stat = statusu)
    if (associated(ta_v)) deallocate(ta_v)
    allocate(ta_v(nv), stat = statusv)
    if (associated(ta_h)) deallocate(ta_h)
    allocate(ta_h(nh), stat = statush)

ta_u=ta_ug(up(:,2)) ! 2: theta-index
!write(*, '(" 2 ")')
!print *, nv, shape(vp), maxval(vp(:,2)), nth, size(ta_vg)
ta_v=ta_vg(vp(:,2))
!write(*, '(" 3 ")')

ta_h=ta_ug(hp(:,2))
!write(*, '(" 4 ")', advance = 'no')

!call CPU_Time(T2)
!     call system_clock ( wall_t2, clock_rate, clock_max )
!     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on u-grid
!write(*, '(" - H,            evaluating on the u-grid ")', advance = 'no')
write(*, '("H, evaluating on the u-, ")', advance = 'no')

!call CPU_Time(T1)
!call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_u)) deallocate(H_u)
    allocate(H_u(nu), stat = statusu)

do cu = 1, nu

      cph=up(cu,1)
      cth=up(cu,2)

      cphr=cph
      cphl=1+modulo((cph-1)-1,nph)

      if (up(cu,3) == 0) then
        H_u(cu) = (H_hg(cphl,cth)+H_hg(cphr,cth))/2
      else
        H_u(cu) = max( H_hg(cphl,cth), H_hg(cphr,cth) )
      endif

enddo

!	SHT-smoothed bottom topography
if ( flag_itd_smooth ) then

	do cu = 1, nu

	      cph=up(cu,1)
	      cth=up(cu,2)

	      cphr=cph
	      cphl=1+modulo((cph-1)-1,nph)

	      if (up(cu,3) == 0) then
	        H_sht_u(cu) = (H_sht_hg(cphl,cth)+H_sht_hg(cphr,cth))/2
	      else
	        H_sht_u(cu) = max( H_sht_hg(cphl,cth), H_sht_hg(cphr,cth) )
	      endif

	enddo

    call save_vector(H_sht_u, dir_cols // 'H_sht_u.dat')
endif
!call CPU_Time(T2)
!     call system_clock ( wall_t2, clock_rate, clock_max )
!     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on v-grid
!write(*, '(" - H,            evaluating on the v-grid ")', advance = 'no')
write(*, '("v-, ")', advance = 'no')

!call CPU_Time(T1)
!call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_v)) deallocate(H_v)
    allocate(H_v(nv), stat = statusv)

do cv = 1, nv

      cph=vp(cv,1)
      cth=vp(cv,2)

      cthb=cth-1
      ctha=cth

      if (vp(cv,3) == 0) then
        H_v(cv) = (H_hg(cph,cthb)+H_hg(cph,ctha))/2
      else
        H_v(cv) = max( H_hg(cph,ctha), H_hg(cph,cthb) )
      endif

enddo

!	SHT-smoothed bottom topography
if ( flag_itd_smooth ) then

	do cv = 1, nv

	      cph=vp(cv,1)
	      cth=vp(cv,2)

	      cthb=cth-1
	      ctha=cth

	      if (vp(cv,3) == 0) then
	        H_sht_v(cv) = (H_sht_hg(cph,cthb)+H_sht_hg(cph,ctha))/2
	      else
	        H_sht_v(cv) = max( H_sht_hg(cph,ctha), H_sht_hg(cph,cthb) )
	      endif

	enddo

    call save_vector(H_sht_v, dir_cols // 'H_sht_v.dat')
endif

!call CPU_Time(T2)
!     call system_clock ( wall_t2, clock_rate, clock_max )
!     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on h-grid
!write(*, '(" - H,            evaluating on the h-grid ")', advance = 'no')
write(*, '("h-grid ")', advance = 'no')

!call CPU_Time(T1)
!call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_h)) deallocate(H_h)
    allocate(H_h(nh), stat = statush)

do ch = 1, nh
      cph=hp(ch,1)
      cth=hp(ch,2)

      H_h(ch) = H_hg(cph, cth)
enddo

!	SHT-smoothed bottom topography
if ( flag_itd_smooth ) then
    call load_alloc_matrix(H_sht_hg, dir_grid // 'H_sht_hg.dat')

	do ch = 1, nh
	      cph=hp(ch,1)
	      cth=hp(ch,2)

	      H_sht_h(ch) = H_sht_hg(cph, cth)
	enddo

    call save_vector(H_sht_h, dir_cols // 'H_sht_h.dat')
endif

!	SHT-smoothed is also used for the barotropic solver
if ( flag_baro_smooth ) then
    call load_alloc_matrix(H_sht_hg, dir_grid // 'H_orig_hg.dat')

	do ch = 1, nh
	      cph=hp(ch,1)
	      cth=hp(ch,2)

	      H_sht_h(ch) = H_sht_hg(cph, cth)
	enddo

    call save_vector(H_sht_h, dir_cols // 'H_orig_h.dat')
endif

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall: '&
						//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')

!%=============================================
!% save columns H_u, H_v, H_h, ta_u, ta_v, ta_h
!%=============================================
     call save_vector(H_u, dir_cols // 'H_u.dat')
     call save_vector(H_v, dir_cols // 'H_v.dat')
     call save_vector(H_h, dir_cols // 'H_h.dat')
     call save_vector(ta_u, dir_cols // 'ta_u.dat')
     call save_vector(ta_v, dir_cols // 'ta_v.dat')
     call save_vector(ta_h, dir_cols // 'ta_h.dat')

end subroutine write_global_ccols_col

!==========================================================================================
! MISC
!==========================================================================================
subroutine metrics_and_masks(nf, P, GD, dir_grid, dir_cols)

!% writes to file the drag of nmodes for all the baro constituents

implicit none

     type(params), intent(in)	:: P
     type(grid_dims), intent(in):: GD

     integer ::  cu,cv,ch,j, nph, nth, nu, nv, nh, nf
     integer :: istat

     character(len = *) :: dir_grid, dir_cols

     real(wp), allocatable		:: metrics_ug(:, :), metrics_vg(:, :)
	 integer, allocatable    	:: up(:, :), vp(:, :), hp(:,:), mask_ug(:,:), mask_vg(:,:), mask_hg(:,:)
	 real(wp), allocatable		:: th_ug(:), ta_ug(:), th_vg(:), ta_vg(:)


! shortcuts
	nph = GD%nph
	nth = GD%nta
	nu = GD%nu
	nv = GD%nv
	nh = GD%nh

	!	nf=ceiling(.5*(P%itm_avg * nph/360 - 1.)) ! will average over 2*nf+1 cells (extra nf in ph and th directions)
	! calculate metrics and masks to later apply running average (u- and v- grids)
	if ( ( (modulo(P%itm_method,10)>=1).or.(P%itd_scheme >= 2) ).and.(nf >= 1) ) then
		if ( (P%messages >= 1).and.(P%smooth_itd_param==1) ) then
			write(*, '("(smoothed over ", a, " grid pts per edge square) ")', advance='no') tostring(2*nf+1)
		endif
	    call load_alloc_matrix(up, dir_cols // 'up.dat')
	    call load_alloc_matrix(vp, dir_cols // 'vp.dat')
	    call load_alloc_matrix(hp, dir_cols // 'hp.dat')

	    ! Mask out the non-ocean grid pts
	    allocate (mask_ug(nph, nth), mask_vg(nph, nth+1), mask_hg(nph, nth), stat=istat)
	    mask_ug = 0; mask_vg = 0; mask_hg = 0;

		do cu = 1, nu
			if (up(cu,3) == 0) then ! mask out boundaries on which u=0
		    	mask_ug(up(cu,1),up(cu,2)) = 1
		    endif
		enddo
		do cv = 1, nv
			if (vp(cv,3) == 0) then ! mask out boundaries on which v=0
		    	mask_vg(vp(cv,1),vp(cv,2)) = 1
		    endif
		enddo
		do ch = 1, nh
			mask_hg(hp(ch,1),hp(ch,2)) = 1
		enddo

	    allocate(metrics_ug(nph, nth), stat = istat)
	    if (P%coor == 1) then
	    	call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
	    	do j = 1,nth
	    		metrics_ug(:, j) = 1 / cosh(ta_ug(j))**2 ! Include metrics for Mercator coords
	    	end do
!	    	deallocate(ta_ug)
	    elseif (P%coor == 2) then
	    	call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
	    	do j = 1,nth
	    		metrics_ug(:, j) = cos(th_ug(j)) ! Include metrics for Lat-Lon coords
	    	end do
!	    	deallocate(th_ug)
	    endif

	    allocate(metrics_vg(nph, nth+1), stat = istat)
	    if (P%coor == 1) then
	    	call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')
	    	do j = 1,nth+1
	    		metrics_vg(:, j) = 1 / cosh(ta_vg(j))**2 ! Include metrics for Mercator coords
	    	end do
	    	deallocate(ta_vg)
	    elseif (P%coor == 2) then
	    	call load_alloc_vector(th_vg, dir_grid // 'th_vg.dat')
	    	do j = 1,nth+1
	    		metrics_vg(:, j) = cos(th_vg(j)) ! Include metrics for Lat-Lon coords
	    	end do
	    	deallocate(th_vg)
	    endif

      call save_matrix(metrics_ug, dir_grid // 'temp/metrics_ug.dat')
      call save_matrix(metrics_vg, dir_grid // 'temp/metrics_vg.dat')
      call save_matrix(mask_ug, dir_grid // 'temp/mask_ug.dat')
      call save_matrix(mask_vg, dir_grid // 'temp/mask_vg.dat')
      call save_matrix(mask_hg, dir_grid // 'temp/mask_hg.dat')

	 endif
	!***************************************************

end subroutine metrics_and_masks

!==========================================================================================

end module generate_matrices




