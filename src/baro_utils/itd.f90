module itd

     use precisions, only: wp, cwp
!     use dispmodule
     use my_trigs
     use my_sparse
     use my_sparse_aggregate
     use control
     use generate_grid, only: grid_dims
     use save_load
     use interpolate
     use my_lapack
     use spline

	 use fd

contains

!**********************************************************************************************************
!**********************************************************************************************************
subroutine eval_c2n_of_H_and_N(P, GD, N_data_dir, hmin, hmax, dir_grid, dir_cols, dir_mats)
!*****************************************************************************
!  eval_c2n_of_H_and_N c2n(H), i.e. squared phase speeds of the normal modes as func of depth
!  for a given a 3D topo Brunt-Vaisala frequency N on a coarse grid (nlon x nlat x nz)
!%================================================================================================
!% 1) Typical atlas grid, like WOA09, is 360x180, then treat each grid cell as a subdomain where
!% 2) Solve the eigenv problem to find c_m(H(x))
!%================================================================================================
!
!     Inputs: 	N = vector of Brunt-Vaisala buoyancy frequency (s^-2)
!		    	z = vector of depth points for N (m)
!				hmin = min depth
!				hmax = max depth
!				xinterp = z coords at which to calc c2n and dc2n
!
!     Outputs:
!                c2n_interp = modal phase speed squared (m/s)
!                dc2n_interp = derivative of the squared modal phase speed (m/s)

implicit none

	type(params), intent(in) :: P
	type(grid_dims), intent(in) :: GD

	real(dp), intent(in)	 :: hmin, hmax

     integer                :: nlon, nlat, nz, clon, clat
     real(wp), allocatable	:: lon(:), lat(:), z(:), dlon_tmp(:), dlat_tmp(:)
     real(wp), allocatable	:: N_topo_3d(:,:,:)
     real(wp)				:: dlon, dlat, H_curr(1), H_curr_orig(1)

	integer					:: nmodes, norder, npl
	real(wp), allocatable	:: cheb_coeff(:,:,:,:) ! two last dims will be allocated as -1,0,1; -1,0,1 with (0,0) corresponding to the current cell and the rest are the neighboring cells
	real(wp), allocatable	:: cheb_coeff_curr(:,:)
!	logical, allocatable	:: interp_mask(:,:)
	real(dp), allocatable	:: H_h(:),H_h_orig(:), H_hg(:,:), H_hg_orig(:,:)
	real(dp), allocatable	:: hmin_local(:,:), hmax_local(:,:)
	integer					:: nu, nv, nh
	real(wp), allocatable	:: c2n_u(:,:), c2n_v(:,:), c2n_h(:,:), dc2n_h(:,:), c2n_h_orig(:,:)
	real(wp), allocatable	:: tmp_u(:), tmp_v(:)
!	real(wp), allocatable	:: wmodes(:,:)
	type (triplet)		    :: im_ph, im_th
	real(wp), allocatable	:: interp_ph(:,:), interp_th(:,:)

	type(triplet)		:: h2u, h2v

	integer					:: nth, nph, cph_neigh, cth_neigh, cph_curr(-1:1), cth_curr(-1:1)
	real(wp), allocatable	:: ph_ug(:), ph_vg(:), th_ug(:), th_vg(:)
!%  upmat(nph, nth): zero if no grid-point, otherwise the u-th point
!%  vpmat(nph, nth+1) : zero if no grid-point, otherwise the v-th point
!%  hpmat(nph, nth) : zero if no grid-point, otherwise the h-th point
	integer, allocatable    :: hpmat(:, :), i_H_h(:)
	real(wp), allocatable   :: sub_c2n_h(:,:), sub_dc2n_h(:,:), sub_c2n_h_orig(:,:)

	integer					:: h,j,k,l, istat!, nz
	integer					:: subd_size_h, sub_nnz_h, cnt_h

	integer					:: iphhl, iphhr, ithhu, ithhl

	integer, allocatable    :: ind_th_ul(:), ind_th_uu(:)
	integer, allocatable    :: ind_ph_vl(:), ind_ph_vr(:)

	character(len = *)		:: dir_cols, dir_grid, N_data_dir, dir_mats

!	real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!	integer :: wall_t1, wall_t2, clock_rate, clock_max

!********************
!	Shortcuts
!********************
	nmodes = P%n_modes
	norder = P%eig_order
!	npl = P%n_cheb

	nph = GD%nph
	nth = GD%nta
	nu = GD%nu
	nv = GD%nv
	nh = GD%nh
!********************
!	Load data
!********************
! Check if WOA09 is being loaded (includes a 3D array with levels)
call load_N_topo_3d(P, N_data_dir, nlon, nlat, nz, lon, lat, z, N_topo_3d)
!N_topo_3d = P%Ns ! test with const stratification

! typically N is on a equally spaced (lat, lon) grid with lon=0.5..359.5; lat=-89.5..89.5
allocate(dlon_tmp(nlon), dlat_tmp(nlat), stat=istat)
dlon_tmp = modulo(cshift(lon, shift=1)-lon, 2*pi)
dlat_tmp = modulo(cshift(lat, shift=1)-lat, pi)
! Check that spacing is uniform. More complicated cases are not yet implemented
if ( ((maxval(dlon_tmp)-minval(dlon_tmp))/minval(dlon_tmp)>1e-6) .or. &
	 ((maxval(dlat_tmp)-minval(dlat_tmp))/minval(dlat_tmp)>1e-6) ) then
	print *, "Horizontal grid for N_topo_3d is non-uniform. This case is not yet implemented, sorry."
	stop
else
	dlon = .5*(maxval(dlon_tmp)+minval(dlon_tmp))
	dlat = .5*(maxval(dlat_tmp)+minval(dlat_tmp))
	if ( (NINT(2*pi/dlon) /= nlon).or.(NINT(pi/dlat) /= nlat) ) then
		print *, "something went wrong..."
	endif
endif

deallocate(dlon_tmp, dlat_tmp)
!********************
!	load grid data
!********************
!	load grid depths
if (P%smooth_type==0) then
	call load_alloc_vector(H_h, dir_cols // 'H_h.dat')
	call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
else
	call load_alloc_vector(H_h, dir_cols // 'H_sht_h.dat')
	call load_alloc_matrix(H_hg, dir_grid // 'H_sht_hg.dat')
	if (P%baro_on_smoothed/=1) then
		call load_alloc_vector(H_h_orig, dir_cols // 'H_h.dat')
		call load_alloc_matrix(H_hg_orig, dir_grid // 'H_hg.dat')
	endif
endif

    call load_alloc_matrix(hpmat, dir_grid // 'ih_hg.dat')

	call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
	call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')

	allocate(ind_th_ul(nlat), ind_th_uu(nlat), stat=istat)
	allocate(ind_ph_vl(nlon), ind_ph_vr(nlon), stat=istat)
!	allocate(ind_ph_u2N(nph), ind_ph_v2N(nph), ind_ph_h2N(nph), stat=istat)
!	allocate(ind_th_u2N(nth), ind_ph_v2N(nth+1), ind_th_h2N(nth), stat=istat)

!	So, N is indeed on a uniform spacial grid. Now match u,v,h-grid points onto corresponding N grid points
ind_ph_vl = 1
ind_ph_vr = nph
do clon = 1,nlon-1
  	call locate(nph, ph_vg, clon*dlon, ind_ph_vr(clon))
  	ind_ph_vl(clon+1) = ind_ph_vr(clon) + 1
enddo

ind_th_ul = 1
ind_th_uu = nth
do clat = 1,nlat
	call locate(nth, th_ug, (-pi/2 + clat*dlat), ind_th_uu(clat))
	if (ind_th_uu(clat) > 0) then
		if (clat < nlat) then
			ind_th_ul(clat+1) = ind_th_uu(clat) + 1
		endif
	else
		ind_th_ul(clat) = 0
		ind_th_uu(clat) = 0
	endif
enddo
where (ind_th_uu<ind_th_ul)
	ind_th_uu=0
	ind_th_ul=0
end where

! Now find min(H_hg) and max(H_hg) on each subdomain. Used in the eigenvalue routine.
	allocate(hmin_local(nlon,nlat), hmax_local(nlon,nlat), stat=istat)
	hmin_local=0.;hmax_local=0.;
do clon = 1,nlon
	do clat = 1,nlat
		!	H  GRID
		iphhl=ind_ph_vl(clon); iphhr=ind_ph_vr(clon); ithhu=ind_th_uu(clat); ithhl=ind_th_ul(clat);
		if ( (iphhl>0).and.(iphhr>0).and.(ithhl>0).and.(ithhu>0) ) then
			if ( maxval(H_hg(iphhl:iphhr, ithhl:ithhu)) >0 ) then
			hmin_local(clon, clat) = minval(H_hg(iphhl:iphhr, ithhl:ithhu), H_hg(iphhl:iphhr, ithhl:ithhu) > 0)
			hmax_local(clon, clat) = maxval(H_hg(iphhl:iphhr, ithhl:ithhu))
		    	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
!				if ( maxval(H_hg_orig(iphhl:iphhr, ithhl:ithhu)) >0 ) then
				hmin_local(clon, clat) = min(hmin_local(clon, clat), minval(H_hg_orig(iphhl:iphhr, ithhl:ithhu), &
																				H_hg_orig(iphhl:iphhr, ithhl:ithhu) > 0) )
				hmax_local(clon, clat) = max(hmax_local(clon, clat), maxval(H_hg_orig(iphhl:iphhr, ithhl:ithhu)) )
!				endif
		    	endif
			endif
		endif

	enddo
enddo
	deallocate(H_hg)
	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
		deallocate(H_hg_orig)
	endif

where (hmin_local > 0) hmin_local = max(hmin_local, hmin)
where (hmax_local > 0) hmax_local = min(hmax_local, hmax)
where (hmax_local <= hmin_local)
	hmax_local = 0
	hmin_local = 0
end where

!*******************************************************
!	Now when subdomains are indexed do the calculations
!*******************************************************

! SWITCH BETWEEN using neighbouring subdomains for bilinear interpolation or not (switch is removed. permanently on)
! Appying bilinear interpolation on precalculated chebychev coeffs results in calc cost increasing by about a factor of 10
! But this should not increase with increasing resolutions. So, probably, not so terrible.
! do bilinear interpolation, nice, smooth c2n but 10 times longer

allocate(c2n_u(nu,nmodes), c2n_v(nv,nmodes), c2n_h(nh,nmodes), dc2n_h(nh,nmodes), stat=istat)
c2n_u = 0; c2n_v = 0; c2n_h = 0; dc2n_h = 0;
if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
	allocate(c2n_h_orig(nh,nmodes), stat=istat)
	c2n_h_orig = 0;
endif

!allocate(interp_mask(-1:1,-1:1), stat=istat)

do clon = 1,nlon
	if (modulo(clon, nlon/10)==0) then
		write(*, '(i3, "%, ")', advance='no') (clon/(nlon/10))*10
	endif

	do clat = 1,nlat
!		print *, clon, clat, hmin_local(clon, clat),hmax_local(clon, clat)
		if ((hmin_local(clon, clat)>0).and.(hmax_local(clon, clat)>0)) then

		! First calc coeffs (for the interval [hmin_local hmax_local])!
		npl = ceiling((hmax_local(clon, clat)-hmin_local(clon, clat))/100) * P%n_cheb ! P%n_cheb points per 1000m of depth

			if (allocated(cheb_coeff)) deallocate(cheb_coeff)
			! two last dims are allocated as -1,0,1; -1,0,1 with (0,0) corresponding to the current cell and the rest are the neighboring cells
			allocate(cheb_coeff(npl,nmodes, -1:1,-1:1), stat=istat)
			cheb_coeff = 0

			if (allocated(cheb_coeff_curr)) deallocate(cheb_coeff_curr)
			allocate(cheb_coeff_curr(npl,nmodes), stat=istat)
			cheb_coeff_curr = 0
!		interp_mask = .false.

!		print *, clon, clat, npl, hmin_local(clon, clat), hmax_local(clon, clat)
			do cph_neigh = -1,1
				cph_curr(cph_neigh) = 1+modulo((clon+cph_neigh)-1,nlon)
				do cth_neigh = -1,1
           			cth_curr(cth_neigh) =  min(max(clat+cth_neigh,1),nlat)
!           			if ((hmin_local(cph_curr(cph_neigh), cth_curr(cth_neigh))>0).and.&
!           				(hmax_local(cph_curr(cph_neigh), cth_curr(cth_neigh))>0)) then
!           				interp_mask(cph_neigh,cth_neigh)=.true.
						!	Warning: this method does not treat the boundary regions any different
						!	Assumes that N_topo_3d is well defined there and can be used for calculations
						call cheb_c2n_vs_depth(hmin_local(clon, clat), hmax_local(clon, clat), &
											 z, N_topo_3d(cph_curr(cph_neigh), cth_curr(cth_neigh),:), nmodes, &
								 				norder, npl, cheb_coeff(:,:, cph_neigh, cth_neigh))
!           			endif
           		enddo
           	enddo

		! Then find grid points at which to eval.

		!	H  GRID
		iphhl=ind_ph_vl(clon); iphhr=ind_ph_vr(clon); ithhu=ind_th_uu(clat); ithhl=ind_th_ul(clat);
		if ( (iphhl>0).and.(iphhr>0).and.(ithhl>0).and.(ithhu>0) ) then

			subd_size_h = (iphhr-iphhl+1)*(ithhu-ithhl+1)
		    sub_nnz_h = count(hpmat(iphhl:iphhr, ithhl:ithhu)>0)

		    if (allocated(i_H_h)) deallocate(i_H_h)
		    allocate(i_H_h(sub_nnz_h), stat = istat)
		    if (allocated(sub_c2n_h)) deallocate(sub_c2n_h)
		    allocate(sub_c2n_h(sub_nnz_h, nmodes), stat = istat)
		    if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
				if (allocated(sub_c2n_h_orig)) deallocate(sub_c2n_h_orig)
		    	allocate(sub_c2n_h_orig(sub_nnz_h, nmodes), stat = istat)
			endif
		    if (allocated(sub_dc2n_h)) deallocate(sub_dc2n_h)
		    allocate(sub_dc2n_h(sub_nnz_h, nmodes), stat = istat)

			!****************************************
	        	call bilinear_2d_mats(3, 3, lon(cph_curr), lat(cth_curr), &
	        				(iphhr-iphhl+1), (ithhu-ithhl+1), ph_vg(iphhl:iphhr), th_ug(ithhl:ithhu), im_ph, im_th)
				call coo2full(im_ph, interp_ph)
				call coo2full(im_th, interp_th)

				cnt_h = 0
			    do k = iphhl,iphhr
					do l = ithhl,ithhu
						if (hpmat(k,l) > 0) then
							cnt_h = cnt_h + 1
							i_H_h(cnt_h) = hpmat(k,l)
							H_curr=H_h(i_H_h(cnt_h))
							if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
								H_curr_orig=H_h_orig(i_H_h(cnt_h))
							endif
!							    	print *, "H_curr = ", H_curr, "H_curr_orig = ", H_curr_orig

							if ( (((P%smooth_type==0).or.(P%smooth_type>=1).and.(P%baro_on_smoothed==1)).and.&
															(H_curr(1) >= hmin).and.(H_curr(1) <= hmax)) .or. &
								(((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)).and.(H_curr(1) >= hmin).and.(H_curr(1) <= hmax).and. &
													(H_curr_orig(1) >= hmin).and.(H_curr_orig(1) <= hmax)) ) then
							  do j = 1, nmodes
								do h = 1, npl
								cheb_coeff_curr(h,j) = vec_vec_dot(interp_th(l-ithhl+1, :), &
												matmul(interp_ph(k-iphhl+1, :),cheb_coeff(h,j, -1:1,-1:1)) )
								enddo

								call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl,&
						                         cheb_coeff_curr(:,j), 1, H_curr, sub_c2n_h(cnt_h,j))
							    if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
									call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl,&
						                         cheb_coeff_curr(:,j), 1, H_curr_orig, sub_c2n_h_orig(cnt_h,j))
								endif
								call chebyshev_deriv_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, &
												cheb_coeff_curr(:,j), 1, H_curr, sub_dc2n_h(cnt_h,j))
						  	  enddo
						  	else
						  	  sub_c2n_h(cnt_h,:) = 0 ! outside of the allowed interval fill with zeros
						  	  if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
						  	  	sub_c2n_h_orig(cnt_h,:) = 0
						  	  endif
						  	  sub_dc2n_h(cnt_h,:) = 0 ! outside of the allowed interval fill with zeros
						  	endif
						endif
					enddo
				enddo

			! now copy the results to the original array
			c2n_h(i_H_h(:), :) = sub_c2n_h(:,:)
			if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
				c2n_h_orig(i_H_h(:), :) = sub_c2n_h_orig(:,:)
			endif
			dc2n_h(i_H_h(:), :) = sub_dc2n_h(:,:)
		endif

		!%===================================
		! Check that everything went OK.
		!%===================================
		!	Potentially c2n and dc2n_h can be negative due to spline extrapolation outside of  [hmin, hmax]. Just chop them off
		where (c2n_h<0) c2n_h = 0
		if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
			where (c2n_h_orig<0) c2n_h_orig = 0
		endif
		where (dc2n_h < 0) dc2n_h = 0

		!	Potentially c2n can be large due to spline extrapolation outside of  [hmin, hmax].
		if ( any( abs(c2n_h(i_H_h(:), 1))>100. ) ) then
			do j = 1, sub_nnz_h
				if ( abs(c2n_h(i_H_h(j), 1))>100) then
					print *, 'Warning!'
					print *, "c2n_h(i_H_h(j), 1),  H_h(i_H_h(j)), hmin_local(clon, clat), hmax_local(clon, clat)"
					print *, c2n_h(i_H_h(j), 1),  H_h(i_H_h(j)), hmin_local(clon, clat), hmax_local(clon, clat)
				endif
			enddo
		endif

		endif ! Refers to: if ((hmin_local(clon, clat)>0).and.(hmax_local(clon, clat)>0)) then
	enddo
enddo

! load sparse matrices
    call load_alloc_sparse(h2u, dir_mats // 'h2u.dat')
	call load_alloc_sparse(h2v, dir_mats // 'h2v.dat')

	allocate(tmp_u(nu), tmp_v(nv), stat=istat)

	do j = 1, nmodes
		call coo_vec_mul(h2u, c2n_h(:,j), P%lib, tmp_u)
		call coo_vec_mul(h2v, c2n_h(:,j), P%lib, tmp_v)
		c2n_u(:,j)=tmp_u
		c2n_v(:,j)=tmp_v
	enddo

	deallocate(tmp_u, tmp_v)
	call deallocate_sparse(h2u)
    call deallocate_sparse(h2v)

	call save_matrix(c2n_u, dir_grid // 'itm/' // 'c2n_u.dat')
	call save_matrix(c2n_v, dir_grid // 'itm/' // 'c2n_v.dat')
	call save_matrix(c2n_h, dir_grid // 'itm/' // 'c2n_h.dat')
	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
		call save_matrix(c2n_h_orig, dir_grid // 'itm/' // 'c2n_h_orig.dat')
	endif
	call save_matrix(dc2n_h, dir_grid // 'itm/' // 'dc2n_h.dat')

end subroutine
!**********************************************************************************************************
!**********************************************************************************************************
subroutine eval_c2n_of_H(P, hmin, hmax, z, N, dir_grid, dir_cols)
!*****************************************************************************
!  eval_c2n_of_H c2n(H), i.e. squared phase speeds of the normal modes as func of depth
!  for a given a column vector of Brunt-Vaisala frequency N as a function of depth z
!  The required modes are found with modecalc_matrix
!
!     Inputs: 	N = vector of Brunt-Vaisala buoyancy frequency (s^-2)
!		    	z = vector of depth points for N (m)
!				hmin = min depth
!				hmax = max depth
!				xinterp = z coords at which to calc c2n and dc2n
!
!     Outputs:
!                c2n_interp = modal phase speed squared (m/s)
!                dc2n_interp = derivative of the squared modal phase speed (m/s)


implicit none

	type(params), intent(in) :: P
	real(dp), intent(in)	:: hmin, hmax
	real(dp), intent(in)	:: N(:), z(:)

	integer					:: nmodes, norder, npl
	real(wp), allocatable	:: cheb_coeff(:,:)
	real(dp), allocatable	:: H_u(:),H_v(:),H_h(:)!, H_h_orig(:)
	integer					:: nu, nv, nh
	real(dp), allocatable	:: c2n_u(:,:), c2n_v(:,:), c2n_h(:,:), dc2n_h(:,:)!, c2n_h_orig(:,:)
!	real(wp), allocatable	:: wmodes(:,:)

	integer j, istat!, nz
	character(len = *)		:: dir_cols, dir_grid

!	real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!	integer :: wall_t1, wall_t2, clock_rate, clock_max


!	shortcuts
	nmodes = P%n_modes

!	load grid depths
if (P%smooth_type==0) then
	call load_alloc_vector(H_u, dir_cols // 'H_u.dat', nu)
	call load_alloc_vector(H_v, dir_cols // 'H_v.dat', nv)
	call load_alloc_vector(H_h, dir_cols // 'H_h.dat', nh)
else
	call load_alloc_vector(H_u, dir_cols // 'H_sht_u.dat', nu)
	call load_alloc_vector(H_v, dir_cols // 'H_sht_v.dat', nv)
	call load_alloc_vector(H_h, dir_cols // 'H_sht_h.dat', nh)

!	call load_alloc_vector(H_h_orig, dir_cols // 'H_h.dat', nh)
endif

	allocate(c2n_u(nu,nmodes), c2n_v(nv,nmodes), c2n_h(nh,nmodes), dc2n_h(nh,nmodes), stat=istat)
!	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
!		allocate(c2n_h_orig(nh,nmodes), stat=istat)
!	endif

if ( size(N)==1) then ! Uniform stratification, use a formula:
	! evaluate c2n on u-, v-grid and dc2n on h-grid
	! c^2_n = (N*H/pi/n)^2
  	do j = 1, nmodes
  		c2n_u(:,j) = (N(1)*H_u/pi/j)**2
  		c2n_v(:,j) = (N(1)*H_v/pi/j)**2
  		c2n_h(:,j) = (N(1)*H_h/pi/j)**2
  		dc2n_h(:,j) = 2*H_h*(N(1)/pi/j)**2
  	enddo

else ! Non-uniform vertical stratification, calc modes numerically
	norder = P%eig_order
	npl = ceiling((hmax-hmin)/100) * P%n_cheb ! P%n_cheb points per 1000m of depth

	allocate(cheb_coeff(npl,nmodes), stat=istat)
	call cheb_c2n_vs_depth(hmin, hmax, z, N, nmodes, norder, npl, cheb_coeff)

	! Interpolate c2n on u-, v-grid and dc2n on h-grid
  do j = 1, nmodes
	call chebyshev_interp ( hmin, hmax, npl, cheb_coeff(:, j), nu, H_u, c2n_u(:, j) )
	call chebyshev_interp ( hmin, hmax, npl, cheb_coeff(:, j), nv, H_v, c2n_v(:, j) )
	call chebyshev_interp ( hmin, hmax, npl, cheb_coeff(:, j), nh, H_h, c2n_h(:, j) )
!	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
!		call chebyshev_interp ( hmin, hmax, npl, cheb_coeff(:, j), nh, H_h_orig, c2n_h_orig(:, j) )
!	endif
	call chebyshev_deriv_interp ( hmin, hmax, npl, cheb_coeff(:, j), nh, H_h, dc2n_h(:, j) )

  	! exclude c2n evaluated at depths less than hmin
	where ((H_u < hmin).or.(H_u > hmax)) c2n_u(:,j) = 0
	where ((H_v < hmin).or.(H_v > hmax)) c2n_v(:,j) = 0
	where ((H_h < hmin).or.(H_h > hmax))
		c2n_h(:,j) = 0
		dc2n_h(:, j) = 0
	end where
  enddo

		!%===================================
		! Check that everything went OK.
		!%===================================
		!	Potentially c2n and dc2n_h can be negative due to spline extrapolation outside of  [hmin, hmax]. Just chop them off
		where (c2n_u<0) c2n_u = 0
		where (c2n_v<0) c2n_v = 0
		where (c2n_h<0) c2n_h = 0
		where (dc2n_h < 0) dc2n_h = 0

endif

	call save_matrix(c2n_u, dir_grid // 'itm/' // 'c2n_u.dat')
	call save_matrix(c2n_v, dir_grid // 'itm/' // 'c2n_v.dat')
	call save_matrix(c2n_h, dir_grid // 'itm/' // 'c2n_h.dat')
	call save_matrix(dc2n_h, dir_grid // 'itm/' // 'dc2n_h.dat')

!	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
!		call save_matrix(c2n_h_orig, dir_grid // 'itm/' // 'c2n_h_orig.dat')
!	endif

end subroutine
!**********************************************************************************************************

subroutine cheb_c2n_vs_depth(hmin, hmax, z, N, nmodes, norder, npl, cheb_coeff)
!*****************************************************************************
!  cheb_c2n_vs_depth c2n(H) calculates the Chebyshev coefficients of the squared phase speeds
!  of the normal modes as func of depth for a given Brunt-Vaisala frequency N as a function of depth z
!  The required modes are found with modecalc_matrix
!
!     Inputs: 	N = vector of Brunt-Vaisala buoyancy frequency (s^-2)
!		    	z = vector of depth points for N (m)
!           	nmodes = number of vertical modes to calculate
!				nord = approx order for modecalc_matrix
!				hmin = min depth
!				hmax = max depth
!				npl = the number of terms in the Chebyshev series
!
!     Outputs:
!                c2n = modal phase speed (m/s)


implicit none

	real(dp), intent(in)	:: hmin, hmax
	real(dp), intent(in)	:: N(:), z(:)
	integer, intent(in)		:: nmodes, norder, npl
	real(wp)				:: xcheb(npl), cheb_coeff(npl,nmodes), c2n(npl,nmodes)
	real(wp), allocatable	:: wmodes(:,:)

	logical, save			:: failed1=.false. ! important! here i rely on keeping the value of failed1 saved between calls of the function
	logical					:: failed2
	real(dp), allocatable	:: N_tmp(:), z_tmp(:)

	integer cz, j, istat, nderiv, npoints, noffset, nz, m, loc

!	real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!	integer :: wall_t1, wall_t2, clock_rate, clock_max


	  if ( size(z) /= size(N) ) then
	    write ( *, '(a)' ) ' '
	    write ( *, '(a)' ) '  Arrays N and z must be of the same size! '
	    stop
	  end if

	! use calc Chebyshev nodes on [hmin, hmax]
	call cheby_zero( hmin, hmax, npl, xcheb )

failed2 = .false.

  	if ( .not.(failed1) ) then
		do j = 1, npl
			call modecalc_matrix(xcheb(j), z, N, nmodes, norder, c2n(j,:), wmodes, failed1) ! careful, wmodes are not calculated
			if (failed1) exit
		enddo
	endif
  	if ( (failed1).and.(.not.(failed2)) ) then
		! Try to double the number of (z, N) points using midpoint interpolation
		call midpoint_interp(z, N, z_tmp, N_tmp, hmin)
!		call disp(z_tmp)
!		call disp(N_tmp)

		do j = 1, npl
			call modecalc_matrix(xcheb(j), z_tmp, N_tmp, nmodes, norder, c2n(j,:), wmodes, failed2)
			if (failed2) then
				print *, "Doubling the number of (z, N) pts didn't help. Stopped."
				stop
			endif
		enddo
	endif

!  do j = 1, npl
!    write ( *, '(2f10.4)' ) xcheb(j), c2n(j,1)
!  end do

  do j = 1, nmodes
	call chebyshev_coefficients( npl, c2n(:,j), cheb_coeff(:, j) )
  enddo

end subroutine

!**********************************************************************************************************

subroutine modecalc_matrix(H, zin, Nin, nmodes, norder, c2n, wmodes, failed)
!*****************************************************************************
!  Modecalc_matrix calculates ocean vertical modes taking
!  a column vector of Brunt-Vaisala frequency N as a function of depth z

!     Inputs:	H = depth for which c2n are calculated
!				N = vector of Brunt-Vaisala buoyancy frequency (s^-2)
!		    	z = vector of depth points (m)
!           	nmodes = number of vertical modes to calculate
!
!     Outputs:   wmodes = vertical velocity structure
!                c2n = squared modal phase speed (m/s)
!
! Eqn for modes: d2w/dz2 + N^2/c^2 *w = 0
! with b.c.: w(-H) = 0 w(0) = 0

implicit none

	real(dp), intent(in)	:: Nin(:), zin(:)
	integer, intent(in)		:: nmodes, norder
	real(wp), intent(in)	:: H
	logical, intent(out)	:: failed

	real(wp)				:: c2n(nmodes)

	real(wp), allocatable	:: N(:), z(:), wmodes(:,:), en(:)

	integer cz, j, istat, nderiv, npoints, noffset, nzin, nz, loc

	real(wp), allocatable	:: fd_coeffs(:,:)
	real(wp), allocatable	:: A(:,:), B(:,:)


!write(*, '("Making the eig. matrix:")')

	nzin = size(zin)
! define the points for the FD scheme
	!	find location of the point deeper than H
		nz = 0
		do j = 1, nzin
	    if (zin(j) - H >= 0) then
	        nz = j
	        exit
	    endif
		end do

	! depth > depth of the last N sample point
		if (nz==0) then
			nz = nzin + CEILING((H - zin(nzin))/(zin(nzin)-zin(nzin-1))) ! extra steps in z
		endif
!	now allocate and fill z, N
	allocate(z(nz),N(nz), stat = istat)
	if (nz<=nzin) then
		z(1:nz-1) = zin(1:nz-1)
		z(nz) = H
		N(1:nz-1) = Nin(1:nz-1)
		! lin interpolate at H
		N(nz) = Nin(nz-1) + (Nin(nz)-Nin(nz-1))/(zin(nz)-zin(nz-1)) * (z(nz)-z(nz-1))
	else
		z(1:nzin) = zin(1:nzin)
		N(1:nzin) = Nin(1:nzin)
		z(nzin+1:nz) = zin(nzin) + (H - zin(nzin))/(nz-nzin) * (/ ( (j), j=1,nz-nzin ) /)
		N(nzin+1:nz) = Nin(nzin) ! fill with the last value of N
	endif
!	print *, z
!	print *, ""
!	print *, N

    ! To approximate Nth deriv with Kth order acc we need N+K+1 points
    nderiv = 2
    npoints = nderiv+norder+1 ! # points to approx the 2nd derivative d2w/dz2
    noffset = floor((npoints-1)/2.)
      if ( npoints > nz ) then
	    write ( *, '(a)' ) ' '
	    print *, nz
	    print *, z(1:nz)
	    write ( *, '(a)' ) '  Not enough Z points to use this FD order! '
	    stop
	  end if


	allocate(fd_coeffs(npoints,nderiv+1), stat = istat)
    allocate(A(nz,nz),B(nz,nz), stat = istat)
    A = 0
    B = 0

    ! calculate the weights for the descrite diff matrix A
    do cz = 2,nz-1
!         select the list of point used
        if (cz <= noffset) then
			call fd_weights (z(cz), z(1:npoints), nderiv, fd_coeffs )
			A(cz, 1:npoints) = fd_coeffs(:, nderiv+1)
        elseif (cz + noffset>=nz) then
			call fd_weights (z(cz), z(nz-npoints+1:nz), nderiv, fd_coeffs )
			A(cz, nz-npoints+1:nz) = fd_coeffs(:, nderiv+1)
        else
			call fd_weights (z(cz), z(cz-noffset:cz-noffset+npoints-1), nderiv, fd_coeffs )
			A(cz, cz-noffset:cz-noffset+npoints-1) = fd_coeffs(:, nderiv+1)
        endif

    enddo
	A = -A ! A ~ -d2/dz2

!   set boundary conditions (zeros at the endpoints)
	A(1,1)=-1.
	A(nz,1)=-1.

!   and B matrix (N^2)
	do cz = 1,nz
		B(cz,cz) = N(cz)**2
	enddo

	allocate(en(nz), stat = istat)
	call gen_eig_nonsym(A,B,nz, en, wmodes) ! wmodes is not calculated unless a corresponding flag is specified

!   extract suitable eigenvalues
!!!!!!!!!!!!!!!! CAREFUL, DOES NOT TAKE CARE OF SELECTING wmodes

!	sort ascending
	call r8vec_sort_bubble_a(nz, en)

!	and find first positive val
	loc = 0
	do j = 1, nz
    if (en(j) >= 1.e-6) then
        loc = j
        exit
    endif
	end do
!print *, "", loc

!	Keep only first nmodes
!	allocate(c2n(nmodes), stat = istat)
		if (loc+nmodes-1 > nz) then
			print *, ""
			print *, "WARNING! Required nmodes is larger than the max number of eigenvals calculated with the given discretization."
			print *, "The number of data points in (z, N) will be doubled (linear midpoint interpolation)."
			failed = .true.
			return
!			stop
		else
			c2n = 1/en(loc:loc+nmodes-1)
		endif
!	print *, cn

end subroutine

!**********************************************************************************************************

subroutine calc_N_z_avrg(P, N_data_dir, dir_grid)

implicit none

	 type(params) :: P
     real(wp), allocatable	:: ph_vg(:), th_ug(:)
     real(wp), allocatable	:: lon(:), lat(:), z(:), H_hg(:,:), H_interp(:,:)
     real(wp), allocatable	:: N_avrg(:), N_topo_3d(:,:,:)
     logical, allocatable	:: mask(:,:)
     real(wp), allocatable	:: dA(:,:)

     integer            :: j, cz, istat, nlon, nlat, nz

     character(len = *)		:: dir_grid, N_data_dir

if (P%N_data .eq. '0') then !no data present
	if (P%N_form == 1) then
		nz = 1
		allocate(z(nz), N_avrg(nz), stat = istat)
		N_avrg = P%Ns
	elseif (P%N_form == 2) then
		nz = 32
		allocate(z(nz), N_avrg(nz), stat = istat)
		z = (/ ( 10*j**2, j=0,nz ) /)
		N_avrg = P%Ns/(1+z/P%NL)
	else
		write(*,'(a)') 'Incorrect choice of input stratification data.'
		stop
	endif
else

	! Check if WOA09 is being loaded (includes a 3D array with levels)
	call load_N_topo_3d(P, N_data_dir, nlon, nlat, nz, lon, lat, z, N_topo_3d)

! load ph_vg, ta_h
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
	 call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')

! interpolate H_hg onto the WOA grid
	call bilinear_2d(size(ph_vg), size(th_ug), ph_vg, th_ug, H_hg, &
    						 nlon, nlat, lon, lat, H_interp, P%lib)

! allocate N_avrg
	allocate(N_avrg(nz), stat = istat)
	allocate(mask(nlon, nlat), stat = istat)
	allocate(dA(nlon, nlat), stat = istat)

		! average horizontally
		do j = 1, nlon
			dA(j, :) = cos(lat)
		enddo

		do cz = 1, nz
			mask = (H_interp > z(cz))
   			N_avrg(cz) = sum(N_topo_3d(:,:,cz)*dA, MASK=mask)/sum(dA, MASK=mask)
    	enddo

	deallocate(dA, mask, H_interp, H_hg, ph_vg, th_ug, N_topo_3d)

endif

call save_vector(z, dir_grid // 'z.dat')
call save_vector(N_avrg, dir_grid // 'N_avrg.dat')


end subroutine calc_N_z_avrg
!--------------------------------------------------------------------------------
subroutine load_N_topo_3d(P, N_data_dir, nlon, nlat, nz, lon, lat, z, N_topo_3d)

implicit none

	 type(params) :: P
     real(wp), allocatable	:: lon(:), lat(:), z(:)
     real(wp), allocatable	:: N_topo_3d(:,:,:)

     integer            :: nlon, nlat, nz

     character(len = *)		:: N_data_dir
     character(len = 100)	:: N_data_file

	! Check if WOA09 is being loaded (includes a 3D array with levels)
	if (scan(P%N_data,'9') > 0) then
		! load N data on the given grid
		N_data_file = N_data_dir // trim(adjustl(P%N_data)) // '.dat'
		call load_alloc_topo_3d(N_data_file, nlon, nlat, nz, lon, lat, z, N_topo_3d, .false.)
		lon = d2r(lon)
		lat = d2r(lat)
	else
		write(*,'(a)') 'WARNING: Please, choose WOA09 stratification file.'
		stop
	endif

end subroutine load_N_topo_3d

!**********************************************************************************************************

subroutine calc_N_on_grid(N_u, N_v, P, nu, nv, nh, N_data_dir, dir_cols, dir_grid, dir_mats)

implicit none

	type(params) :: P
     real(wp), allocatable	:: ph_vg(:), th_ug(:)
     real(wp), allocatable	:: lon(:), lat(:), z(:), H_hg(:,:)
     integer, allocatable	:: hp(:,:)
     real(wp), allocatable	:: N_topo(:,:), N_hg(:,:)
     real(wp), allocatable	:: N_topo_3d(:,:,:), N_hg_3d(:,:,:)
     type (triplet)			:: h2v, h2u

     integer, intent(in)    :: nu, nv, nh
     integer            :: j, cz, istat, nlon, nlat, nz

     character(len = *)		:: dir_grid, dir_cols, dir_mats, N_data_dir
     character(len = 100)	:: N_data_file

     real(wp), allocatable, dimension(:) :: N_u, N_v, N_h

! allocate N_u and N_v
    allocate(N_u(nu),N_v(nv), N_h(nh), stat = istat)

if (P%N_data .eq. '0') then !no data present
	N_u = P%Ns
	N_v = P%Ns
else
! load hp, ph_vg, ta_h
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')

! load sparse matrices
    call load_alloc_sparse(h2u, dir_mats // 'h2u.dat')
	call load_alloc_sparse(h2v, dir_mats // 'h2v.dat')

! load N data on the given grid
	N_data_file = N_data_dir // trim(adjustl(P%N_data)) // '.dat'
	! Check if WOA09 is being loaded (includes a 3D array with levels)
	if (scan(P%N_data,'9') > 0) then
		call load_alloc_topo_3d(N_data_file, nlon, nlat, nz, lon, lat, z, N_topo_3d, .false.)
		allocate(N_hg_3d(size(ph_vg), size(th_ug), nz), stat = istat)
		! interpolate horizontally onto our grid
		do cz = 1, nz
    		call bilinear_2d(nlon, nlat, d2r(lon), d2r(lat), N_topo_3d(:,:,cz), &
    						 size(ph_vg), size(th_ug), ph_vg, th_ug, N_topo, P%lib)
			N_hg_3d(:,:,cz) = N_topo
    	enddo
    	! now interpolate vertically
		if (P%smooth_type==0) then
			call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
		else
			call load_alloc_matrix(H_hg, dir_grid // 'H_sht_hg.dat')
		endif
    	call interp_z_3d_to_2d(size(ph_vg), size(th_ug), nz, z, H_hg, N_hg_3d, N_hg)
    	deallocate(N_topo, N_topo_3d, N_hg_3d, H_hg)
	else
	! Assume WOA05 is being loaded (vertically averaged, 2D array)
		call load_alloc_topo(N_data_file, nlon, nlat, lon, lat, N_topo, .false.)
		! interpolate on our grid
    	call bilinear_2d(nlon, nlat, d2r(lon), d2r(lat), N_topo, size(ph_vg), size(th_ug), ph_vg, th_ug, N_hg, P%lib)
    	deallocate(N_topo)
	endif

! Because of spline interpolation there is a few negative values of N_hg
! Just turn them positive. So few that it doesn't matter
	N_hg = abs(N_hg)
!call save_matrix(N_hg, dir_grid // 'N_hg.dat')

 ! write N_u, N_v, N_h in columns
     N_u = 0
     N_v = 0
     N_h = 0

    do j = 1, nh
		N_h(j) = N_hg(hp(j,1),hp(j,2))
	end do
	call coo_vec_mul(h2u, N_h, P%lib, N_u)
	call coo_vec_mul(h2v, N_h, P%lib, N_v)

    call deallocate_sparse(h2u)
    call deallocate_sparse(h2v)
    deallocate(hp, ph_vg, th_ug, N_h)

end if

end subroutine calc_N_on_grid
!**********************************************************************************************************
subroutine interp_z_3d_to_2d(nlon, nlat, nz, z, H_hg, N_topo_3d, N_topo)

	implicit none

	integer	:: nlon, nlat, nz
	real ( kind = 8 ), allocatable	:: ypp(:)
	real ( kind = 8 )				:: ypval, yppval, zero_dp = 0

	real(wp), intent(in)	:: z(nz), N_topo_3d(:,:,:), H_hg(:,:)
	real(wp), allocatable	:: N_topo(:,:)
	real ( kind = 8 )		:: N_tmp(nlon,nlat)
	integer :: istat, clon, clat


    allocate(ypp(nz), stat = istat)

do clon = 1,nlon
	do clat = 1,nlat
!	 set the spline to have 0 derivative at the right boundary
!						and 0 second derivative at the left endpoint
		call spline_cubic_set ( nz, real(z,kind=8), real(N_topo_3d(clon,clat,:),kind=8), 2, zero_dp, 1, zero_dp, ypp ) !0 - the spline should be a quadratic over the first and last intervals
		call spline_cubic_val ( nz, real(z,kind=8), real(N_topo_3d(clon,clat,:),kind=8), ypp, &
														real(H_hg(clon,clat),kind=8), N_tmp(clon,clat), ypval, yppval )
	enddo
enddo

! back to the variables in our precision
allocate(N_topo(nlon,nlat), stat = istat)
!overwrite extrapolated values for the deep ocean with deepest given N
where (H_hg > z(nz))
	N_topo = N_topo_3d(:,:,nz)
elsewhere
	N_topo = N_tmp
end where

end subroutine interp_z_3d_to_2d

!**********************************************************************************************************
subroutine midpoint_interp(z, N, z_interp, N_interp, zlim)
!

	implicit none

	real(dp), intent(in)	:: z(:), N(:), zlim
	real(dp), allocatable	:: z_interp(:), N_interp(:)

	integer :: istat, c_interp, niz, nz, j, loc

		nz = size(z)
	!	find location of the point deeper than H
		loc = 0
		do j = 1, nz
		    if (z(j) - zlim >= 0) then
		        loc = j
		        exit
		    endif
		end do
	! depth > depth of the last N sample point
		if (loc==0) then
			loc = nz
		endif

!	niz = 2*(nz-1)+1 ! don't double the # of pts
	niz = 2*(loc-1)+1 + (nz-loc) ! don't double the # of pts between z(1) and z(loc)
	if (.not.(allocated(z_interp))) allocate(z_interp(niz),N_interp(niz), stat = istat)
    z_interp = 0
    N_interp = 0

	! copy vals for the original pts
	do c_interp = 1,loc
		z_interp(1 + 2*(c_interp-1)) = z(c_interp)
		N_interp(1 + 2*(c_interp-1)) = N(c_interp)
	enddo
	if (nz>loc) then
		z_interp(2*loc:niz) = z(loc+1:nz) ! the outside pts remain the same
		N_interp(2*loc:niz) = N(loc+1:nz)
	endif

	! midpoint interpolate vals inbetween
	do c_interp = 2,2*(loc-1),2
		z_interp(c_interp) = (z_interp(c_interp+1) + z_interp(c_interp-1))/2
		N_interp(c_interp) = (N_interp(c_interp+1) + N_interp(c_interp-1))/2
	enddo

end subroutine midpoint_interp
!**********************************************************************************************************

subroutine calc_itd_prefactor(ufactor,vfactor, N_u, N_v, dHdph, dHdta, cpt, P, nu, nv, dir_cols, dir_grid)

implicit none

	type(params) :: P
     real(wp), allocatable	:: N_u(:), N_v(:),  dHdph(:), dHdta(:), igw_slope_u(:), igw_slope_v(:)
     real(wp), allocatable	:: ta_vg(:), ph_vg(:), ta_ug(:), ph_ug(:)
     real(wp), allocatable	:: ta_v(:), ph_v(:), ta_u(:), ph_u(:), th_v(:), th_u(:)
     real(wp), allocatable	:: ph_u_nonrot(:), th_u_nonrot(:), ph_v_nonrot(:), th_v_nonrot(:)
	 real(wp)           :: th0, ph0
	 real(wp)           :: area_ufactor, area_vfactor
	 integer, allocatable   :: up(:,:), vp(:,:)

     real(wp), allocatable	:: issponge_u(:, :), issponge_v(:, :)

     integer, intent(in)    :: nu, nv!, nh
     integer            :: j, istat

	 character(len=2)   :: cpt
	 type (tide_params) :: pars
     character(len = *) :: dir_grid, dir_cols!, dir_mats, N_data_dir

     integer, allocatable, dimension(:) :: ufactor,vfactor


allocate(ufactor(nu), vfactor(nv), stat = istat)
ufactor = 0
vfactor = 0

if (P%trapped == 0) then !no do not cut off ITD conversion above crit lats
	ufactor = 1
	vfactor = 1
else
! load hp, ph_vg, ta_h
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')

     call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')

    allocate(ta_v(nv), ph_v(nv), ta_u(nu), ph_u(nu), th_v(nv), th_u(nu), stat = istat)
! and fill them

	ph_u = (/ ( (ph_ug(up(j, 1))), j=1,nu ) /)
	ph_v = (/ ( (ph_vg(vp(j, 1))), j=1,nv ) /)
	ta_u = (/ ( (ta_ug(up(j, 2))), j=1,nu ) /)
	ta_v = (/ ( (ta_vg(vp(j, 2))), j=1,nv ) /)

	th_u=ta2th(ta_u, P%coor)
	th_v=ta2th(ta_v, P%coor)

	deallocate(ph_ug, ph_vg, ta_ug, ta_vg, ta_u, ta_v)
! poles
th0 = d2r(P%latP)
ph0 = d2r(P%lonP)
! transform to the non-rotated coords
	call calc_latlon_nonrot(ph_u, th_u, ph0, th0, ph_u_nonrot, th_u_nonrot)
	call calc_latlon_nonrot(ph_v, th_v, ph0, th0, ph_v_nonrot, th_v_nonrot)

	pars = get_pars(cpt)

	where ((abs(th_u_nonrot)) < asin(.5*pars%omega0/P%omega)) ufactor = 1
	where ((abs(th_v_nonrot)) < asin(.5*pars%omega0/P%omega))  vfactor = 1

	if (P%itd_scheme == 3) then

		!  keep parametrized ITD drag only outside of the mode 1 sponge region
		if (P%sponge_itd_param == 0) then
			call load_alloc_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//cpt(2:2)//'.dat')
			call load_alloc_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//cpt(2:2)//'.dat')

			where (issponge_u(:,1)/=0)
				ufactor = 0
			end where
			where (issponge_v(:,1)/=0)
				vfactor = 0
			end where
		endif


	! keep parametrized ITD drag only at "rough" topography (where  eps=|grad(H)|/tan(alpha) > eps_rough )
		area_ufactor = dot_product(ufactor,cos(th_u))
		area_vfactor = dot_product(vfactor,cos(th_v))
		! calculate eps=|grad(H)|/tan(alpha)
		allocate(igw_slope_u(nu), igw_slope_v(nv), stat = istat)

		igw_slope_u = real( ( ((pars%omega0)**2 - (2*P%omega*sin(th_u_nonrot))**2)/(N_u**2-(pars%omega0)**2) )**0.5 )
		igw_slope_v = real( ( ((pars%omega0)**2 - (2*P%omega*sin(th_v_nonrot))**2)/(N_v**2-(pars%omega0)**2) )**0.5 )

		where (abs(dHdph) < P%eps_rough * igw_slope_u) ufactor = 0 ! topography is not rough enough
		where (abs(dHdta) < P%eps_rough * igw_slope_v) vfactor = 0 ! topography is not rough enough

		if (P%messages == 2) then
			print *, ""
			write(*, '("param ITD for ", a, " is used over (% area): ", f5.1, ", ", f5.1, " ")') cpt, &
				100*dot_product(ufactor,cos(th_u))/area_ufactor, 100*dot_product(vfactor,cos(th_v))/area_vfactor
		endif
	endif


	call save_vector(ufactor, dir_cols // 'ufactor_'//cpt//'.dat')
	call save_vector(vfactor, dir_cols // 'vfactor_'//cpt//'.dat')

!	print *, sum(ufactor), size(ufactor)
!	print *, sum(vfactor), size(vfactor)

	deallocate(ph_u, th_u,ph_v, th_v, ph_u_nonrot, th_u_nonrot, ph_v_nonrot, th_v_nonrot)

end if

end subroutine calc_itd_prefactor

!**********************************************************************************************************
subroutine calc_crit_lat_prefactor(u_crit_lat,v_crit_lat, cpts, P, nu, nv, dir_cols, dir_grid)

implicit none

	type(params) :: P
     real(wp), allocatable	:: ta_vg(:), ph_vg(:), ta_ug(:), ph_ug(:)
     real(wp), allocatable	:: ta_v(:), ph_v(:), ta_u(:), ph_u(:), th_v(:), th_u(:)
     real(wp), allocatable	:: ph_u_nonrot(:), th_u_nonrot(:), ph_v_nonrot(:), th_v_nonrot(:)
	 real(wp)           :: th0, ph0

	 integer, allocatable   :: up(:,:), vp(:,:)

     integer, intent(in)    :: nu, nv
     integer            :: j, istat, ncpts, ccpt

     character(len=*) :: cpts
	 character(len=2)   :: cpt
	 type (tide_params) :: pars
     character(len = *) :: dir_grid, dir_cols!, dir_mats, N_data_dir

     logical, allocatable, dimension(:,:) :: u_crit_lat,v_crit_lat


! load hp, ph_vg, ta_h
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')

     call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')

    allocate(ta_v(nv), ph_v(nv), ta_u(nu), ph_u(nu), th_v(nv), th_u(nu), stat = istat)
! and fill them

	ph_u = (/ ( (ph_ug(up(j, 1))), j=1,nu ) /)
	ph_v = (/ ( (ph_vg(vp(j, 1))), j=1,nv ) /)
	ta_u = (/ ( (ta_ug(up(j, 2))), j=1,nu ) /)
	ta_v = (/ ( (ta_vg(vp(j, 2))), j=1,nv ) /)

	th_u=ta2th(ta_u, P%coor)
	th_v=ta2th(ta_v, P%coor)

	deallocate(up, vp, ph_ug, ph_vg, ta_ug, ta_vg, ta_u, ta_v)

	! poles
	th0 = d2r(P%latP)
	ph0 = d2r(P%lonP)
	! transform to the non-rotated coords
	call calc_latlon_nonrot(ph_u, th_u, ph0, th0, ph_u_nonrot, th_u_nonrot)
	call calc_latlon_nonrot(ph_v, th_v, ph0, th0, ph_v_nonrot, th_v_nonrot)


ncpts=len(cpts)/2

allocate(u_crit_lat(ncpts,nu), v_crit_lat(ncpts,nv), stat = istat)
! .true. within the (-lat_critical, lat_critical) where ITs radiate away
u_crit_lat = .false.
v_crit_lat = .false.

	do ccpt = 1, ncpts

	    cpt=cpts(2*ccpt-1:2*ccpt)

		pars = get_pars(cpt)

		where ((abs(th_u_nonrot)) < asin(.5*pars%omega0/P%omega)) u_crit_lat(ccpt, :) = .true.
		where ((abs(th_v_nonrot)) < asin(.5*pars%omega0/P%omega))  v_crit_lat(ccpt, :) = .true.

	enddo

!	call save_vector(u_crit_lat, dir_cols // 'u_crit_lat_'//cpt//'.dat')
!	call save_vector(v_crit_lat, dir_cols // 'v_crit_lat_'//cpt//'.dat')

	deallocate(ph_u, th_u,ph_v, th_v, ph_u_nonrot, th_u_nonrot, ph_v_nonrot, th_v_nonrot)

end subroutine calc_crit_lat_prefactor


!**********************************************************************************************************
subroutine eval_c2n_of_H_and_N_old(P, N_data_dir, hmin, hmax, dir_grid, dir_cols)
!*****************************************************************************
!  eval_c2n_of_H_and_N c2n(H), i.e. squared phase speeds of the normal modes as func of depth
!  for a given a 3D topo Brunt-Vaisala frequency N on a coarse grid (nlon x nlat x nz)
!%================================================================================================
!% 1) Typical atlas grid, like WOA09, is 360x180, then treat each grid cell as a subdomain where
!% 2) Solve the eigenv problem to find c_m(H(x))
!%================================================================================================
!
!     Inputs: 	N = vector of Brunt-Vaisala buoyancy frequency (s^-2)
!		    	z = vector of depth points for N (m)
!				hmin = min depth
!				hmax = max depth
!				xinterp = z coords at which to calc c2n and dc2n
!
!     Outputs:
!                c2n_interp = modal phase speed squared (m/s)
!                dc2n_interp = derivative of the squared modal phase speed (m/s)

implicit none

	type(params), intent(in) :: P
	real(dp), intent(in)	 :: hmin, hmax

     integer                :: nlon, nlat, nz, clon, clat
     real(wp), allocatable	:: lon(:), lat(:), z(:), dlon_tmp(:), dlat_tmp(:)
     real(wp), allocatable	:: N_topo_3d(:,:,:)
     real(wp)				:: dlon, dlat, H_curr(1), H_curr_orig(1)

	integer					:: nmodes, norder, npl
	real(wp), allocatable	:: cheb_coeff(:,:,:,:) ! two last dims will be allocated as -1,0,1; -1,0,1 with (0,0) corresponding to the current cell and the rest are the neighboring cells
	real(wp), allocatable	:: cheb_coeff_curr(:,:)
!	logical, allocatable	:: interp_mask(:,:)
	real(dp), allocatable	:: H_u(:),H_v(:),H_h(:),H_h_orig(:), H_hg(:,:), H_ug(:,:), H_vg(:,:), H_hg_orig(:,:)
	real(dp), allocatable	:: hmin_local(:,:), hmax_local(:,:)
	integer					:: nu, nv, nh, cu, cv
	real(wp), allocatable	:: c2n_u(:,:), c2n_v(:,:), c2n_h(:,:), dc2n_h(:,:), c2n_h_orig(:,:)
!	real(wp), allocatable	:: wmodes(:,:)
	type (triplet)		    :: im_ph, im_th
	real(wp), allocatable	:: interp_ph(:,:), interp_th(:,:)
	integer					:: switch ! switch on/off calcs on neighbours for bilinear interpolation

	integer					:: nth, nph, cph_neigh, cth_neigh, cph_curr(-1:1), cth_curr(-1:1)
	real(wp), allocatable	:: ph_ug(:), ph_vg(:), th_ug(:), th_vg(:)
!%  upmat(nph, nth): zero if no grid-point, otherwise the u-th point
!%  vpmat(nph, nth+1) : zero if no grid-point, otherwise the v-th point
!%  hpmat(nph, nth) : zero if no grid-point, otherwise the h-th point
	integer, allocatable    :: upmat(:, :), vpmat(:, :), hpmat(:, :)
	integer, allocatable    :: up(:, :), vp(:, :)
	integer, allocatable    :: i_H_u(:), i_H_v(:), i_H_h(:)
	real(wp), allocatable   :: sub_c2n_u(:,:), sub_c2n_v(:,:), sub_c2n_h(:,:), sub_dc2n_h(:,:), sub_c2n_h_orig(:,:)

	integer					:: h,j,k,l, istat!, nz
	integer					:: subd_size_u, subd_size_v, subd_size_h
	integer					:: sub_nnz_u, sub_nnz_v, sub_nnz_h
	integer					:: cnt_u,cnt_v,cnt_h

	integer					:: iphul, iphur, ithuu, ithul
	integer					:: iphvl, iphvr, ithvu, ithvl
	integer					:: iphhl, iphhr, ithhu, ithhl

	integer, allocatable    :: ind_ph_ul(:), ind_ph_ur(:), ind_th_ul(:), ind_th_uu(:)
	integer, allocatable    :: ind_ph_vl(:), ind_ph_vr(:), ind_th_vl(:), ind_th_vu(:)

	character(len = *)		:: dir_cols, dir_grid, N_data_dir

!	real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!	integer :: wall_t1, wall_t2, clock_rate, clock_max

!********************
!	Shortcuts
!********************
	nmodes = P%n_modes
	norder = P%eig_order
!	npl = P%n_cheb
!********************
!	Load data
!********************
! Check if WOA09 is being loaded (includes a 3D array with levels)
call load_N_topo_3d(P, N_data_dir, nlon, nlat, nz, lon, lat, z, N_topo_3d)
!N_topo_3d = P%Ns ! test with const stratification

! typically N is on a equally spaced (lat, lon) grid with lon=0.5..359.5; lat=-89.5..89.5
allocate(dlon_tmp(nlon), dlat_tmp(nlat), stat=istat)
dlon_tmp = modulo(cshift(lon, shift=1)-lon, 2*pi)
dlat_tmp = modulo(cshift(lat, shift=1)-lat, pi)
! Check that spacing is uniform. More complicated cases are not yet implemented
if ( ((maxval(dlon_tmp)-minval(dlon_tmp))/minval(dlon_tmp)>1e-6) .or. &
	 ((maxval(dlat_tmp)-minval(dlat_tmp))/minval(dlat_tmp)>1e-6) ) then
	print *, "Horizontal grid for N_topo_3d is non-uniform. This case is not yet implemented, sorry."
	stop
else
	dlon = .5*(maxval(dlon_tmp)+minval(dlon_tmp))
	dlat = .5*(maxval(dlat_tmp)+minval(dlat_tmp))
	if ( (NINT(2*pi/dlon) /= nlon).or.(NINT(pi/dlat) /= nlat) ) then
		print *, "something went wrong..."
	endif
endif

deallocate(dlon_tmp, dlat_tmp)
!********************
!	load grid data
!********************
!	load grid depths
if (P%smooth_type==0) then
	call load_alloc_vector(H_u, dir_cols // 'H_u.dat', nu)
	call load_alloc_vector(H_v, dir_cols // 'H_v.dat', nv)
	call load_alloc_vector(H_h, dir_cols // 'H_h.dat', nh)
	call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
else
	call load_alloc_vector(H_u, dir_cols // 'H_sht_u.dat', nu)
	call load_alloc_vector(H_v, dir_cols // 'H_sht_v.dat', nv)
	call load_alloc_vector(H_h, dir_cols // 'H_sht_h.dat', nh)
	call load_alloc_matrix(H_hg, dir_grid // 'H_sht_hg.dat')
	if (P%baro_on_smoothed/=1) then
		call load_alloc_vector(H_h_orig, dir_cols // 'H_h.dat', nh)
		call load_alloc_matrix(H_hg_orig, dir_grid // 'H_hg.dat')
	endif
endif

     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(upmat, dir_grid // 'iu_ug.dat')
     call load_alloc_matrix(vpmat, dir_grid // 'iv_vg.dat')
     call load_alloc_matrix(hpmat, dir_grid // 'ih_hg.dat')

	call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat', nph)
	call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
	call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat', nth)
	call load_alloc_vector(th_vg, dir_grid // 'th_vg.dat')

	allocate(H_ug(nph,nth),H_vg(nph,nth+1), stat=istat)
	H_ug = 0; H_vg = 0;
	do cu = 1, nu
	    H_ug(up(cu,1),up(cu,2)) = H_u(cu)
	enddo
	do cv = 1, nv
	    H_vg(vp(cv,1),vp(cv,2)) = H_v(cv)
	enddo
	deallocate(up, vp)

	allocate(ind_ph_ul(nlon), ind_ph_ur(nlon), ind_th_ul(nlat), ind_th_uu(nlat), stat=istat)
	allocate(ind_ph_vl(nlon), ind_ph_vr(nlon), ind_th_vl(nlat), ind_th_vu(nlat), stat=istat)
!	allocate(ind_ph_u2N(nph), ind_ph_v2N(nph), ind_ph_h2N(nph), stat=istat)
!	allocate(ind_th_u2N(nth), ind_ph_v2N(nth+1), ind_th_h2N(nth), stat=istat)

!	So, N is indeed on a uniform spacial grid. Now match u,v,h-grid points onto corresponding N grid points
!	First for u-grid
ind_ph_ul = 1
ind_ph_ur = nph
do clon = 1,nlon-1
  	call locate(nph, ph_ug, clon*dlon, ind_ph_ur(clon))
  	ind_ph_ul(clon+1) = ind_ph_ur(clon) + 1
enddo
ind_th_ul = 1
ind_th_uu = nth
do clat = 1,nlat
	call locate(nth, th_ug, (-pi/2 + clat*dlat), ind_th_uu(clat))
	if (ind_th_uu(clat) > 0) then
		if (clat < nlat) then
			ind_th_ul(clat+1) = ind_th_uu(clat) + 1
		endif
	else
		ind_th_ul(clat) = 0
		ind_th_uu(clat) = 0
	endif
enddo
where (ind_th_uu<ind_th_ul)
	ind_th_uu=0
	ind_th_ul=0
end where
!	Then for v-grid
ind_ph_vl = 1
ind_ph_vr = nph
do clon = 1,nlon-1
  	call locate(nph, ph_vg, clon*dlon, ind_ph_vr(clon))
  	ind_ph_vl(clon+1) = ind_ph_vr(clon) + 1
enddo
ind_th_vl = 1
ind_th_vu = size(th_vg)
do clat = 1,nlat
	call locate(nth, th_vg, (-pi/2 + clat*dlat), ind_th_vu(clat))
	if (ind_th_vu(clat) > 0) then
		if (clat < nlat) then
			ind_th_vl(clat+1) = ind_th_vu(clat) + 1
		endif
	else
		ind_th_vl(clat) = 0
		ind_th_vu(clat) = 0
	endif
enddo
where (ind_th_vu<ind_th_vl)
	ind_th_vu=0
	ind_th_vl=0
end where

! Now find min(H_hg) and max(H_hg) on each subdomain. Used in the eigenvalue routine.
	allocate(hmin_local(nlon,nlat), hmax_local(nlon,nlat), stat=istat)
	hmin_local=0.
	hmax_local=0.
do clon = 1,nlon
	do clat = 1,nlat
		!	U  GRID
		iphul=ind_ph_ul(clon); iphur=ind_ph_ur(clon); ithuu=ind_th_uu(clat); ithul=ind_th_ul(clat);
		if ( (iphul>0).and.(iphur>0).and.(ithul>0).and.(ithuu>0) ) then
			if ( maxval(H_ug(iphul:iphur, ithul:ithuu)) >0 ) then
				hmin_local(clon, clat) = minval(H_ug(iphul:iphur, ithul:ithuu), H_ug(iphul:iphur, ithul:ithuu) > 0)
				hmax_local(clon, clat) = maxval(H_ug(iphul:iphur, ithul:ithuu))
			endif
		endif
		!	V  GRID
		iphvl=ind_ph_vl(clon); iphvr=ind_ph_vr(clon); ithvu=ind_th_vu(clat); ithvl=ind_th_vl(clat);
		if ( (iphvl>0).and.(iphvr>0).and.(ithvl>0).and.(ithvu>0) ) then
			if ( maxval(H_vg(iphvl:iphvr, ithvl:ithvu)) >0 ) then
			hmin_local(clon, clat) = min(hmin_local(clon, clat), minval(H_vg(iphvl:iphvr, ithvl:ithvu), H_vg(iphvl:iphvr, ithvl:ithvu) > 0) )
			hmax_local(clon, clat) = max(hmax_local(clon, clat), maxval(H_vg(iphvl:iphvr, ithvl:ithvu)) )
			endif
		endif
		!	H  GRID
		iphhl=ind_ph_vl(clon); iphhr=ind_ph_vr(clon); ithhu=ind_th_uu(clat); ithhl=ind_th_ul(clat);
		if ( (iphhl>0).and.(iphhr>0).and.(ithhl>0).and.(ithhu>0) ) then
			if ( maxval(H_hg(iphhl:iphhr, ithhl:ithhu)) >0 ) then
			hmin_local(clon, clat) = min(hmin_local(clon, clat), minval(H_hg(iphhl:iphhr, ithhl:ithhu), H_hg(iphhl:iphhr, ithhl:ithhu) > 0) )
			hmax_local(clon, clat) = max(hmax_local(clon, clat), maxval(H_hg(iphhl:iphhr, ithhl:ithhu)) )
		    	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
!				if ( maxval(H_hg_orig(iphhl:iphhr, ithhl:ithhu)) >0 ) then
				hmin_local(clon, clat) = min(hmin_local(clon, clat), minval(H_hg_orig(iphhl:iphhr, ithhl:ithhu), &
																				H_hg_orig(iphhl:iphhr, ithhl:ithhu) > 0) )
				hmax_local(clon, clat) = max(hmax_local(clon, clat), maxval(H_hg_orig(iphhl:iphhr, ithhl:ithhu)) )
!				endif
		    	endif
			endif
		endif

	enddo
enddo
	deallocate(H_ug, H_vg, H_hg)
	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
		deallocate(H_hg_orig)
	endif

where (hmin_local > 0) hmin_local = max(hmin_local, hmin)
where (hmax_local > 0) hmax_local = min(hmax_local, hmax)
where (hmax_local <= hmin_local)
	hmax_local = 0
	hmin_local = 0
end where

!call save_matrix(hmin_local, dir_grid // 'hmin_local.dat')
!call save_matrix(hmax_local, dir_grid // 'hmax_local.dat')

!	print *, ""
!	print *, ind_ph_vl(350:360)
!	print *, ind_ph_vr(350:360)
!	print *, ind_th_vl(170:180)
!	print *, ind_th_vu(170:180)

!*******************************************************
!	Now when subdomains are indexed do the calculations
!*******************************************************

! SWITCH BETWEEN using neighbouring subdomains for bilinear interpolation or not
! Appying bilinear interpolation on precalculated chebychev coeffs results in calc cost increasing by about a factor of 10
! But this should not increase with increasing resolutions. So, probably, not so terrible.
!switch = 0 ! original version. no interpolation. jumps in cn at the boundarier between subdomains
switch = 1 ! do bilinear interpolation, nice but 10 times longer

allocate(c2n_u(nu,nmodes), c2n_v(nv,nmodes), c2n_h(nh,nmodes), dc2n_h(nh,nmodes), stat=istat)
c2n_u = 0; c2n_v = 0; c2n_h = 0; dc2n_h = 0;
if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
	allocate(c2n_h_orig(nh,nmodes), stat=istat)
	c2n_h_orig = 0;
endif

!allocate(interp_mask(-1:1,-1:1), stat=istat)

do clon = 1,nlon
	if (modulo(clon, nlon/10)==0) then
		write(*, '(i3, "%, ")', advance='no') (clon/(nlon/10))*10
	endif

	do clat = 1,nlat
!		print *, clon, clat
		if ((hmin_local(clon, clat)>0).and.(hmax_local(clon, clat)>0)) then

		! First calc coeffs (for the interval [hmin_local hmax_local])!
		npl = ceiling((hmax_local(clon, clat)-hmin_local(clon, clat))/100) * P%n_cheb ! P%n_cheb points per 1000m of depth

		if (switch == 1) then
			if (allocated(cheb_coeff)) deallocate(cheb_coeff)
			! two last dims are allocated as -1,0,1; -1,0,1 with (0,0) corresponding to the current cell and the rest are the neighboring cells
			allocate(cheb_coeff(npl,nmodes, -1:1,-1:1), stat=istat)
			cheb_coeff = 0

			if (allocated(cheb_coeff_curr)) deallocate(cheb_coeff_curr)
			allocate(cheb_coeff_curr(npl,nmodes), stat=istat)
			cheb_coeff_curr = 0
!		interp_mask = .false.

!		print *, clon, clat, npl, hmin_local(clon, clat), hmax_local(clon, clat)
			do cph_neigh = -1,1
				cph_curr(cph_neigh) = 1+modulo((clon+cph_neigh)-1,nlon)
				do cth_neigh = -1,1
           			cth_curr(cth_neigh) =  min(max(clat+cth_neigh,1),nlat)
!           			if ((hmin_local(cph_curr(cph_neigh), cth_curr(cth_neigh))>0).and.&
!           				(hmax_local(cph_curr(cph_neigh), cth_curr(cth_neigh))>0)) then
!           				interp_mask(cph_neigh,cth_neigh)=.true.
						!	Warning: this method does not treat the boundary regions any different
						!	Assumes that N_topo_3d is well defined there and can be used for calculations
						call cheb_c2n_vs_depth(hmin_local(clon, clat), hmax_local(clon, clat), &
											 z, N_topo_3d(cph_curr(cph_neigh), cth_curr(cth_neigh),:), nmodes, &
								 				norder, npl, cheb_coeff(:,:, cph_neigh, cth_neigh))
!           			endif
           		enddo
           	enddo
        else
			if (allocated(cheb_coeff_curr)) deallocate(cheb_coeff_curr)
			allocate(cheb_coeff_curr(npl,nmodes), stat=istat)
			cheb_coeff_curr = 0

			call cheb_c2n_vs_depth(hmin_local(clon, clat), hmax_local(clon, clat), z, N_topo_3d(clon, clat, :),&
						  			nmodes, norder, npl, cheb_coeff_curr(:,:))
        endif

		! Then find grid points at which to eval.
		!	U  GRID
		iphul=ind_ph_ul(clon); iphur=ind_ph_ur(clon); ithuu=ind_th_uu(clat); ithul=ind_th_ul(clat);

		if ( (iphul>0).and.(iphur>0).and.(ithul>0).and.(ithuu>0) ) then

			subd_size_u = (iphur-iphul+1)*(ithuu-ithul+1)
		    sub_nnz_u = count(upmat(iphul:iphur, ithul:ithuu)>0)

		    if (allocated(i_H_u)) deallocate(i_H_u)
		    allocate(i_H_u(sub_nnz_u), stat = istat)
		    if (allocated(sub_c2n_u)) deallocate(sub_c2n_u)
		    allocate(sub_c2n_u(sub_nnz_u, nmodes), stat = istat)

			if (switch == 1) then
			!****************************************
	        	call bilinear_2d_mats(3, 3, lon(cph_curr), lat(cth_curr), &
	        				(iphur-iphul+1), (ithuu-ithul+1), ph_ug(iphul:iphur), th_ug(ithul:ithuu), im_ph, im_th)
				call coo2full(im_ph, interp_ph)
				call coo2full(im_th, interp_th)

!			if (clon==nlon) then
!				call disp(interp_ph)
!				call disp(interp_th)
!				stop
!			endif

				cnt_u = 0
			    do k = iphul,iphur
					do l = ithul,ithuu
						if (upmat(k,l) > 0) then
							cnt_u = cnt_u + 1
							i_H_u(cnt_u) = upmat(k,l)
	!						write ( *, '(3f10.3)' ) interp_th(l-ithul+1, :)
	!						write ( *, '(3f10.3)' ) interp_ph(k-iphul+1, :)
							H_curr=H_u(i_H_u(cnt_u))

							if ((H_curr(1) >= hmin).and.(H_curr(1) <= hmax)) then
							  do j = 1, nmodes
								do h = 1, npl
								cheb_coeff_curr(h,j) = vec_vec_dot(interp_th(l-ithul+1, :), &
												matmul(interp_ph(k-iphul+1, :),cheb_coeff(h,j, -1:1,-1:1)) )
!								call disp(cheb_coeff(h,j, -1:1,-1:1))
!								call disp(cheb_coeff_curr(h,j))
								enddo

								call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl,&
						                         cheb_coeff_curr(:,j), 1, H_curr, sub_c2n_u(cnt_u,j))
						  	  enddo
						  	else
						  	  sub_c2n_u(cnt_u,:) = 0 ! outside of the allowed interval fill with zeros
						  	endif
						endif
					enddo
				enddo
			!****************************************
			else
			!****************************************
				cnt_u = 0
			    do k = iphul,iphur
					do l = ithul,ithuu
						if (upmat(k,l) > 0) then
							cnt_u = cnt_u + 1
							i_H_u(cnt_u) = upmat(k,l)
						endif
					enddo
				enddo
			! Interpolate c2n on u-grid
			  do j = 1, nmodes
				  call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, cheb_coeff_curr(:,j),&
			                        	sub_nnz_u, H_u(i_H_u(:)), sub_c2n_u(:,j))
				  ! exclude c2n evaluated at depths less than hmin
				  where ((H_u(i_H_u(:)) < hmin).or.(H_u(i_H_u(:)) > hmax)) sub_c2n_u(:,j) = 0
			  enddo
			!****************************************
			endif
		! now copy the results to the original array
		  c2n_u(i_H_u(:), :) = sub_c2n_u(:,:)

		endif

		!	V  GRID
		iphvl=ind_ph_vl(clon); iphvr=ind_ph_vr(clon); ithvu=ind_th_vu(clat); ithvl=ind_th_vl(clat);
		if ( (iphvl>0).and.(iphvr>0).and.(ithvl>0).and.(ithvu>0) ) then

			subd_size_v = (iphvr-iphvl+1)*(ithvu-ithvl+1)
		    sub_nnz_v = count(vpmat(iphvl:iphvr, ithvl:ithvu)>0)

		    if (allocated(i_H_v)) deallocate(i_H_v)
		    allocate(i_H_v(sub_nnz_v), stat = istat)
		    if (allocated(sub_c2n_v)) deallocate(sub_c2n_v)
		    allocate(sub_c2n_v(sub_nnz_v, nmodes), stat = istat)

			if (switch == 1) then
			!****************************************
	        	call bilinear_2d_mats(3, 3, lon(cph_curr), lat(cth_curr), &
	        				(iphvr-iphvl+1), (ithvu-ithvl+1), ph_vg(iphvl:iphvr), th_vg(ithvl:ithvu), im_ph, im_th)
				call coo2full(im_ph, interp_ph)
				call coo2full(im_th, interp_th)

				cnt_v = 0
			    do k = iphvl,iphvr
					do l = ithvl,ithvu
						if (vpmat(k,l) > 0) then
							cnt_v = cnt_v + 1
							i_H_v(cnt_v) = vpmat(k,l)
							H_curr=H_v(i_H_v(cnt_v))

							if ((H_curr(1) >= hmin).and.(H_curr(1) <= hmax)) then
							  do j = 1, nmodes
								do h = 1, npl
								cheb_coeff_curr(h,j) = vec_vec_dot(interp_th(l-ithvl+1, :), &
												matmul(interp_ph(k-iphvl+1, :),cheb_coeff(h,j, -1:1,-1:1)) )
								enddo

								call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl,&
						                         cheb_coeff_curr(:,j), 1, H_curr, sub_c2n_v(cnt_v,j))
						  	  enddo
						  	else
						  	  sub_c2n_v(cnt_v,:) = 0 ! outside of the allowed interval fill with zeros
						  	endif
						endif
					enddo
				enddo
			!****************************************
			else
			!****************************************
				cnt_v = 0
			    do k = iphvl,iphvr
					do l = ithvl,ithvu
						if (vpmat(k,l) > 0) then
							cnt_v = cnt_v + 1
							i_H_v(cnt_v) = vpmat(k,l)
						endif
					enddo
				enddo

			! Interpolate c2n on u-grid
			  do j = 1, nmodes
				call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, cheb_coeff_curr(:,j), &
				                        sub_nnz_v, H_v(i_H_v(:)), sub_c2n_v(:,j))
				  ! exclude c2n evaluated at depths less than hmin
				  where ((H_v(i_H_v(:)) < hmin).or.(H_v(i_H_v(:)) > hmax)) sub_c2n_v(:,j) = 0
			  enddo
			!****************************************
			endif

			! now copy the results to the original array
			c2n_v(i_H_v(:), :) = sub_c2n_v(:,:)
		endif

		!	H  GRID
		iphhl=ind_ph_vl(clon); iphhr=ind_ph_vr(clon); ithhu=ind_th_uu(clat); ithhl=ind_th_ul(clat);
		if ( (iphhl>0).and.(iphhr>0).and.(ithhl>0).and.(ithhu>0) ) then

			subd_size_h = (iphhr-iphhl+1)*(ithhu-ithhl+1)
		    sub_nnz_h = count(hpmat(iphhl:iphhr, ithhl:ithhu)>0)

		    if (allocated(i_H_h)) deallocate(i_H_h)
		    allocate(i_H_h(sub_nnz_h), stat = istat)
		    if (allocated(sub_c2n_h)) deallocate(sub_c2n_h)
		    allocate(sub_c2n_h(sub_nnz_h, nmodes), stat = istat)
		    if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
				if (allocated(sub_c2n_h_orig)) deallocate(sub_c2n_h_orig)
		    	allocate(sub_c2n_h_orig(sub_nnz_h, nmodes), stat = istat)
			endif
		    if (allocated(sub_dc2n_h)) deallocate(sub_dc2n_h)
		    allocate(sub_dc2n_h(sub_nnz_h, nmodes), stat = istat)

			if (switch == 1) then
			!****************************************
	        	call bilinear_2d_mats(3, 3, lon(cph_curr), lat(cth_curr), &
	        				(iphhr-iphhl+1), (ithhu-ithhl+1), ph_vg(iphhl:iphhr), th_ug(ithhl:ithhu), im_ph, im_th)
				call coo2full(im_ph, interp_ph)
				call coo2full(im_th, interp_th)

				cnt_h = 0
			    do k = iphhl,iphhr
					do l = ithhl,ithhu
						if (hpmat(k,l) > 0) then
							cnt_h = cnt_h + 1
							i_H_h(cnt_h) = hpmat(k,l)
							H_curr=H_h(i_H_h(cnt_h))
							if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
								H_curr_orig=H_h_orig(i_H_h(cnt_h))
							endif
!							    	print *, "H_curr = ", H_curr, "H_curr_orig = ", H_curr_orig

							if ( (((P%smooth_type==0).or.(P%smooth_type>=1).and.(P%baro_on_smoothed==1)).and.&
															(H_curr(1) >= hmin).and.(H_curr(1) <= hmax)) .or. &
								(((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)).and.(H_curr(1) >= hmin).and.(H_curr(1) <= hmax).and. &
													(H_curr_orig(1) >= hmin).and.(H_curr_orig(1) <= hmax)) ) then
							  do j = 1, nmodes
								do h = 1, npl
								cheb_coeff_curr(h,j) = vec_vec_dot(interp_th(l-ithhl+1, :), &
												matmul(interp_ph(k-iphhl+1, :),cheb_coeff(h,j, -1:1,-1:1)) )
								enddo

								call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl,&
						                         cheb_coeff_curr(:,j), 1, H_curr, sub_c2n_h(cnt_h,j))
							    if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
									call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl,&
						                         cheb_coeff_curr(:,j), 1, H_curr_orig, sub_c2n_h_orig(cnt_h,j))
								endif
								call chebyshev_deriv_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, &
												cheb_coeff_curr(:,j), 1, H_curr, sub_dc2n_h(cnt_h,j))
						  	  enddo
						  	else
						  	  sub_c2n_h(cnt_h,:) = 0 ! outside of the allowed interval fill with zeros
						  	  if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
						  	  	sub_c2n_h_orig(cnt_h,:) = 0
						  	  endif
						  	  sub_dc2n_h(cnt_h,:) = 0 ! outside of the allowed interval fill with zeros
						  	endif
						endif
					enddo
				enddo
			!****************************************
			else
			!****************************************
				cnt_h = 0
			    do k = iphhl,iphhr
					do l = ithhl,ithhu
						if (hpmat(k,l) > 0) then
							cnt_h = cnt_h + 1
							i_H_h(cnt_h) = hpmat(k,l)
						endif
					enddo
				enddo

				! Interpolate c2n and dc2n on h-grid
				  do j = 1, nmodes
					call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, cheb_coeff_curr(:,j), &
					                        sub_nnz_h, H_h(i_H_h(:)), sub_c2n_h(:,j))
					if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
						call chebyshev_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, cheb_coeff_curr(:,j), &
						                    sub_nnz_h,  H_h_orig(i_H_h(:)), sub_c2n_h_orig(:,j))
					endif
					call chebyshev_deriv_interp ( hmin_local(clon, clat), hmax_local(clon, clat), npl, cheb_coeff_curr(:,j), &
					                              sub_nnz_h, H_h(i_H_h(:)), sub_dc2n_h(:, j) )

					  ! exclude c2n evaluated at depths less than hmin
					    if ( (P%smooth_type==0).or.(P%smooth_type>=1).and.(P%baro_on_smoothed==1) ) then
						  where ( (H_h(i_H_h(:)) < hmin).or.(H_h(i_H_h(:)) > hmax) )
						    sub_c2n_h(:,j) = 0
							sub_dc2n_h(:, j) = 0
					  	  end where
					    else
						  where ( ((H_h(i_H_h(:)) < hmin).or.(H_h(i_H_h(:)) > hmax)) .or. &
					    		  ((H_h_orig(i_H_h(:)) < hmin).or.(H_h_orig(i_H_h(:)) > hmax)) )
						    sub_c2n_h(:,j) = 0
							sub_dc2n_h(:, j) = 0
					  	  end where
					    endif
				  enddo
			  !****************************************
			  endif

			! now copy the results to the original array
			c2n_h(i_H_h(:), :) = sub_c2n_h(:,:)
			if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
				c2n_h_orig(i_H_h(:), :) = sub_c2n_h_orig(:,:)
			endif
			dc2n_h(i_H_h(:), :) = sub_dc2n_h(:,:)
		endif

		!%===================================
		! Check that everything went OK.
		!%===================================
		!	Potentially c2n and dc2n_h can be negative due to spline extrapolation outside of  [hmin, hmax]. Just chop them off
		where (c2n_u<0) c2n_u = 0
		where (c2n_v<0) c2n_v = 0
		where (c2n_h<0) c2n_h = 0
		if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
			where (c2n_h_orig<0) c2n_h_orig = 0
		endif
		where (dc2n_h < 0) dc2n_h = 0

		!	Potentially c2n can be large due to spline extrapolation outside of  [hmin, hmax].
		if (any( abs(c2n_u(i_H_u(:), 1))>100. ) .or. any( abs(c2n_v(i_H_v(:), 1))>100. )) then
			do j = 1, sub_nnz_u
				if ( abs(c2n_u(i_H_u(j), 1))>100) then
					print *, 'Warning!'
					print *, "c2n_u(i_H_u(j), 1),  H_u(i_H_u(j)), hmin_local(clon, clat), hmax_local(clon, clat)"
					print *, c2n_u(i_H_u(j), 1),  H_u(i_H_u(j)), hmin_local(clon, clat), hmax_local(clon, clat)
				endif
			enddo
			do j = 1, sub_nnz_v
				if ( abs(c2n_v(i_H_v(j), 1))>100) then
					print *, 'Warning!'
					print *, "c2n_v(i_H_v(j), 1),  H_v(i_H_u(j)), hmin_local(clon, clat), hmax_local(clon, clat)"
					print *, c2n_v(i_H_v(j), 1),  H_v(i_H_v(j)), hmin_local(clon, clat), hmax_local(clon, clat)
				endif
			enddo
		endif

		endif ! Refers to: if ((hmin_local(clon, clat)>0).and.(hmax_local(clon, clat)>0)) then
	enddo
enddo

	call save_matrix(c2n_u, dir_grid // 'itm/' // 'c2n_u.dat')
	call save_matrix(c2n_v, dir_grid // 'itm/' // 'c2n_v.dat')
	call save_matrix(c2n_h, dir_grid // 'itm/' // 'c2n_h.dat')
	if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
		call save_matrix(c2n_h_orig, dir_grid // 'itm/' // 'c2n_h_orig.dat')
	endif
	call save_matrix(dc2n_h, dir_grid // 'itm/' // 'dc2n_h.dat')

end subroutine
!**********************************************************************************************************


end module itd
