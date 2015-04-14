module generate_grid

     use precisions, only: dp, eps
     use my_trigs
     use err_manager
     use dispmodule
	 use control
     use save_load
	 use read_etopo
	 use interpolate
	 use arrays
	 use fd
     use spherepack_iface

     use testshrotate

     implicit none

!********************************************************************

type grid_dims

     integer ::     coor ! 1 for Mercator, 2 for lat-lon
     integer ::     nph, nta
     integer ::     nu=0, nv=0, nh=0, np=0
     real(dp)::		lonP, latP

end type

     contains
!***********************************************************************************************************
subroutine generate_global_grid(etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols, P, GD)

     character(len=*) :: etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols
     type(params)	  :: P
     type(grid_dims)  :: GD

     integer :: numLons, numLats
     real(dp), allocatable :: xValues(:), yValues(:)
     integer, allocatable  :: zValues(:, :)
     real(dp), allocatable :: ph(:), th(:)!, topo(:)

     real(dp), allocatable :: ph_vg(:), th_ug(:) !ta_ug(:), ta_vg(:),
     real(dp), allocatable :: H_hg(:, :), H_sht_hg(:,:)

     real(dp), allocatable	:: smooth(:)
     integer				:: j,k
     real(dp), allocatable	:: H_sht(:,:), topo_sht(:,:)
     real(dp), allocatable	:: lon_sht(:), lat_sht(:)
     integer			    :: nlon_sht, nlat_sht, istat

     character(len=10)	:: str_tmp
     character(len=50)	:: filename
     character(len=100) :: dirname
     logical			:: load_smooth_file

     real    ::      T1, T2! for measuring CPU (NOT REAL TIME!)
     integer :: 	 wall_t1, wall_t2, clock_rate, clock_max

	real(dp), allocatable :: ph_tmp(:), th_tmp(:)
	integer					:: nph_tmp, nth_tmp
!%================================================
!% Loading/preparing topography file
!%================================================

     if (P%load_etopo==1) then
		numLons = P%etopo_res;
		numLats = P%etopo_res/2;
		call prepare_topo(etopo_file, numLons, numLats, xValues, yValues, zValues, &
                            real(P%lonP, kind=dp), real(P%latP, kind=dp))
		! add to the filename the resolution in minutes
		write (str_tmp, '(g10.2)') 360*60/real(numLons-1)
		filename = 'topo_rot_'//trim(adjustl(str_tmp))//'min_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat'
		dirname = topo_dir_out
		!write to the file
		call save_topo(topo_dir_out // trim(filename), numLons, numLats, xValues, yValues, zValues)

     else
		  filename = topo_file
		  dirname = topo_dir_in
          call load_alloc_topo(topo_dir_in//topo_file, numLons, numLats, xValues, yValues, zValues, (P%messages>=1))
     endif

! convert xValues and yValues to radians
    allocate(ph(size(xValues)), th(size(yValues)), stat = istat)
    ph=d2r(xValues)
    th=d2r(yValues)
	deallocate(xValues, yValues)

!%===================================================
!*********** IRRELEVANT TO THE MODEL ***********
!% Apply S-G filter directly. For spectrum comparison
!%===================================================
if (1==0) then
	nph_tmp = P%nph+1
	nth_tmp = P%nph/2+1
    allocate(ph_tmp(nph_tmp))
    allocate(th_tmp(nth_tmp))

	ph_tmp=(/ ( (2*pi*j/P%nph), j=0,nph_tmp-1 ) /)
	th_tmp=(/ ( (-pi/2 + 2*pi*j/P%nph), j=0,nth_tmp-1 ) /)

!% Interpolate bathymetry to global grid H_hg
    call bilinear_2d(numLons, numLats, ph, th, real(zValues,dp), nph_tmp, nth_tmp, ph_tmp, th_tmp, H_hg, P%lib)

	! save interpolated bathy
	write (str_tmp, '(g10.2)') 360*60/real(nph_tmp-1)
	filename = 'topo_rot_'//trim(adjustl(str_tmp))//'min_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'_interp.dat'
	call save_topo(topo_dir_in // trim(filename), nph_tmp, nth_tmp, r2d(ph_tmp), r2d(th_tmp), nint(H_hg))

	allocate(H_sht_hg(nph_tmp,nth_tmp), stat=istat)
	! I want the topography to vary quadratically within two min wavelengths
	! which is why I use pts: -P%igw_gppw:P%igw_gppw in x and y to obtain a quadratic approx for the topo
	call savgol_smooth_2d(nph_tmp,nth_tmp, H_hg, P%igw_gppw, P%igw_gppw, 2, 2, H_sht_hg)

	! save SG filtered bathy
	filename = 'topo_rot_'//trim(adjustl(str_tmp))//'min_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'_sg.dat'
	call save_topo(topo_dir_in // trim(filename), nph_tmp, nth_tmp, r2d(ph_tmp), r2d(th_tmp), nint(H_sht_hg))

	stop
endif
!%================================================
!% set up the global C-grid and operator matrices
!%================================================
	call gen_global_grid(zValues, ph, th, numLons-1, H_hg, dir_grid, P)
      ! (ph_ug, ph_vg, ta_ug, ta_vg) are saved into /global/cols

  GD%coor = P%coor
  GD%lonP = P%lonP
  GD%latP = P%latP
  GD%nph = P%nph
  GD%nta = size(H_hg,2)

      if ((minval(H_hg(:,1)) < 0.).or.(minval(H_hg(:,size(H_hg,2))) < 0.)) then
	      write(*, '("****** Ocean is in the first/last tau index of H_hg. Removed ******")')
    	    where (H_hg(:,1)<0) H_hg(:,1)=real(0,dp)! new depth is 0
	        where (H_hg(:,size(H_hg,2))<0) H_hg(:,size(H_hg,2))=real(0,dp)! new depth is 0
	  end if

      ! set the minimum permitted water depth on the global grid:
      call adjust_min_depth(H_hg, GD%nph, GD%nta, P, dir_grid)
!call save_matrix(H_hg, dir_grid // 'H_hg_raw.dat')

      ! bays
      if (P%bays == 1) then
      	call remove_bays(H_hg, GD%nph, GD%nta, real(0,dp))!, dir_grid)  ! new depth is 0
      end if
!call save_matrix(H_hg, dir_grid // 'H_hg_raw1.dat')
      ! seas
      if (P%inlandseas == 1) then
		write(*, '("Removing inland seas/lakes: ")', advance = 'no')
      	call remove_seas(H_hg, GD%nph, GD%nta, real(0,dp), GD%nph/10)  ! new depth is 0
      end if
!call save_matrix(H_hg, dir_grid // 'H_hg_raw2.dat')
!     Make H_hg POSITIVE, until here below the sea level H_hg was NEGATIVE.
      H_hg = -H_hg
!     save H on the C-grid
      call save_matrix(H_hg, dir_grid // 'H_hg.dat')

!%============================================
! Smooth bottom topography for the ITD scheme
!%============================================
if ( (P%itd_scheme>0).and.(P%smooth_type>0).and.(P%sht_smooth_H >= 0) ) then
	call CPU_Time(T1)
	call system_clock ( wall_t1, clock_rate, clock_max )

	write(*, '("Smoothing topo for the linear ITD scheme. Smoothing type - ")', advance = 'no')

	!%==========================
	if (P%smooth_type<=3) then ! Use a global topo filter, based on SH transform
	!%==========================
		! IMPORTANT: given value of sht_smooth_H is used only with itd_scheme=1.
		! When itd_scheme>1, smoothing degree is calculated using igw_gppw parameter
		if (P%itd_scheme>=2) then
			! largest grip spacing is: dx = 2*pi*re/nph; smallest SH wavelength: lambda = pi*re/n
			! lambda = gppw*dx  => n = nph/(2*gppw)
			if (P%smooth_type==1) then! 1 - truncation smoothing
				P%sht_smooth_H =P%nph/(2*P%igw_gppw)
				write(*, '("SH truncation; ")', advance = 'no')
			elseif (P%smooth_type==2) then ! 2 - Gaussian smoothing
				P%sht_smooth_H =P%nph/(2*P%igw_gppw)
				write(*, '("Gaussian smoothing wrt SH; ")', advance = 'no')
			elseif (P%smooth_type==3) then
				P%sht_smooth_H =P%nph/(P%igw_gppw) ! 3 - truncation + sinc smoothing; (increse the # of harmonics by two since the second half is heavily averaged out)
				write(*, '("sinc + truncation filter wrt SH ; ")', advance = 'no')
			endif
		endif
		write(*, '("n_sht=", a,"; ")', advance = 'no') tostring(P%sht_smooth_H)

!		Attempt to load the smoothed topography

		inquire( file=trim(dirname)//'smoothed_t'//tostring(P%smooth_type)//'_trunc'//tostring(P%sht_smooth_H)// '_'// &
   	               				 trim(filename), exist=load_smooth_file )
		if (load_smooth_file) then
   			call load_alloc_topo(trim(dirname)//'smoothed_t'//tostring(P%smooth_type)//'_trunc'//tostring(P%sht_smooth_H)// '_'// &
   	               				 trim(filename), nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht, (P%messages>=1))
		else
	     	call sht_smooth_H(ph, th,real(zValues,dp), P, lon_sht, lat_sht, H_sht, trim(dirname), trim(filename))
	    endif
		! DEALLOCATE TOPO
		deallocate(ph, th, zValues)

     	! shortcuts
		nlon_sht = size(lon_sht)
		nlat_sht = size(lat_sht)
		call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
		call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
	! Now restrict the smoothed topography (for ITD) onto the grid topo H_hg (interpolate)
	    call bilinear_2d(nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht, &
	                                         GD%nph, GD%nta, ph_vg, th_ug, H_sht_hg, P%lib)
		deallocate(lon_sht, lat_sht, H_sht, ph_vg, th_ug)

	  if (P%baro_on_smoothed==1) then
	      ! set the minimum permitted water depth on the global grid:
	      call adjust_min_depth(H_sht_hg, GD%nph, GD%nta, P, dir_grid)
	  endif

!     Make H_sht_hg POSITIVE, until here below the sea level H_hg was NEGATIVE.
      	H_sht_hg = -H_sht_hg

	!%==========================
	elseif (P%smooth_type==4) then ! Apply Savitzky-Galoy filter
	!%==========================
		write(*, '("Savitzky-Galoy filter; n_pts=", a, ", n_ord=", i1,"; ")', advance = 'no') tostring(P%igw_gppw), 2
		allocate(H_sht_hg(GD%nph, GD%nta), stat=istat)
		! I want the topography to vary quadratically within two min wavelengths
		! which is why I use pts: -P%igw_gppw:P%igw_gppw in x and y to obtain a quadratic approx for the topo
		call savgol_smooth_2d(GD%nph, GD%nta, H_hg, P%igw_gppw, P%igw_gppw, 2, 2, H_sht_hg)

		if (P%baro_on_smoothed==1) then
			H_sht_hg = -H_sht_hg
			! set the minimum permitted water depth on the global grid:
	      	call adjust_min_depth(H_sht_hg, GD%nph, GD%nta, P, dir_grid)
	        H_sht_hg = -H_sht_hg
	    endif
	endif
	!%==========================

		! Ad-hoc: fix continent boundaries
		where ( (H_hg <= 0) ) H_sht_hg = 0
		! Ad-hoc: do not want another grid indexing for smoothed topography. So, for points where (H_sht_hg<=0) but (H_hg > 0) use vals of H_hg
		! typically along the continental boundaries (i.e., the Gibbs phenomenon)
		where ( (H_hg > 0).and.(H_sht_hg<=0) ) H_sht_hg = H_hg

        call save_matrix(H_sht_hg, dir_grid // 'H_sht_hg.dat')
        deallocate(H_sht_hg)

	call CPU_Time(T2)
	call system_clock ( wall_t2, clock_rate, clock_max )
	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall: '&
						//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')
endif

	!%================================================================================================
	!% If the barotropic problem is to be solved on a smoothed grid, then link H_hg to H_sht_hg
	!%================================================================================================
	if ((P%smooth_type>0).and.(P%baro_on_smoothed)==1) then
		call system('mv ' // dir_grid//'H_hg.dat ' // dir_grid//'H_orig_hg.dat')
		call system('ln -s ' // dir_grid//'H_sht_hg.dat ' // dir_grid//'H_hg.dat')
	endif

!     print *,  P%nph, GD%nta, P%nph*GD%nta, count(H_hg>0), count(H_hg>0)/(P%nph*GD%nta)
     call write_global_cgrid(H_hg, GD%nph, GD%nta, GD%nu, GD%nv, GD%nh, GD%np, dir_grid, dir_cols)

end subroutine

!***********************************************************************************************************

subroutine write_global_cgrid(H, nph, nth, nu, nv, nh, np, dir_grid, dir_cols)

    implicit none
!% for a height field H(nph, nth), defined on the h-grid,
!% generates appropriate u, v and h grids in an Arakawa C manner.
!% Returns :
!%
!%  up(1:nu,1:4) : phi-index, theta-index, boundary-flag, opt. position
!%  vp(1:nv,1:4) : phi-index, theta-index, boundary-flag, opt. position
!%  hp(1:nh,1:4) : phi-index, theta-index,             0, opt. position
!%
!%  upmat(nph, nth): zero if no grid-point, otherwise the u-th point
!%  vpmat(nph, nth+1) : zero if no grid-point, otherwise the v-th point
!%  hpmat(nph, nth) : zero if no grid-point, otherwise the h-th point

     integer, intent(in)     ::     nph, nth
     real(dp), allocatable:: H(:, :)
     integer, allocatable :: up(:, :), vp(:, :), hp(:, :)

     integer, allocatable :: upmat(:, :), vpmat(:, :), hpmat(:, :)
     integer ::     nu, nv, nh, np

     integer ::      statusu, statusv, statush
     integer ::      cth, cph, cphl, cphr

     character(len=*) :: dir_grid, dir_cols

write(*, '("Allocating grid points: ")', advance = 'no')

    if (allocated(up)) deallocate(up)
    allocate(up(nph*nth,4), stat = statusu)
    if (allocated(vp)) deallocate(vp)
    allocate(vp(nph*nth,4), stat = statusv)
    if (allocated(hp)) deallocate(hp)
    allocate(hp(nph*nth,4), stat = statush)

up = 0
vp = 0
hp = 0

    if (allocated(upmat)) deallocate(upmat)
    allocate(upmat(nph,nth), stat = statusu)
    if (allocated(vpmat)) deallocate(vpmat)
    allocate(vpmat(nph,nth+1), stat = statusv)
    if (allocated(hpmat)) deallocate(hpmat)
    allocate(hpmat(nph,nth), stat = statush)

upmat = 0
vpmat = 0
hpmat = 0

do cph = 1,nph
!print *, ""
!write(*, '(i3)') cph
     do cth = 1,nth
!write(*, '(i3,  " ")', advance = 'no') cth
           cphl=1+modulo((cph-1)-1,nph)
           cphr=1+modulo((cph+1)-1,nph)

         if (H(cph, cth) > 0) then
         if ( (H(cphr, cth) > 0) .or. (H(cphl, cth) > 0) .or. &
                (H(cph, cth-1) > 0) .or. (H(cph, cth+1) > 0) ) then

             !% allot an h-grid point:
             nh=nh+1
             np=np+1

             hp(nh,1)=cph
             hp(nh,2)=cth
             hp(nh,4)=np
             hpmat(cph, cth)=nh

             !% ... and some u-grid points:
             nu=nu+1
             np=np+1
             up(nu,1)=cph
             up(nu,2)=cth
             up(nu,4)=np
             upmat(cph, cth)=nu

             if (H(cphl, cth) <= 0) then
               up(nu,3)=1
             endif

             if (H(cphr, cth) <= 0) then
               nu=nu+1
               np=np+1
               up(nu,1)=cphr
               up(nu,2)=cth
               up(nu,3)=1
               up(nu,4)=np
               upmat(cphr, cth)=nu
             endif

             !% .. and some v-grid points:
             nv=nv+1
             np=np+1
             vp(nv,1)=cph
             vp(nv,2)=cth
             vp(nv,4)=np
             vpmat(cph,cth)=nv
             if (H(cph, cth-1) <= 0) then
               vp(nv,3)=1
             endif
             if (H(cph, cth+1) <= 0) then
               nv=nv+1
               np=np+1
               vp(nv,1)=cph
               vp(nv,2)=cth+1
               vp(nv,3)=1
               vp(nv,4)=np
               vpmat(cph,cth+1)=nv
             endif

         endif
         endif

     enddo
enddo

write(*, '( a, " for u, ", a, " for v, ", a, "  for h.")') tostring(nu),tostring(nv),tostring(nh)

!%=============================================
!% save the forward and backward C grid indexes
!%=============================================
! save iu_ug, iv_vg, ih_hg
     call save_matrix(upmat, dir_grid // 'iu_ug.dat')
     call save_matrix(vpmat, dir_grid // 'iv_vg.dat')
     call save_matrix(hpmat, dir_grid // 'ih_hg.dat')
! save up, vp, hp
     call save_matrix(up(1:nu,1:4), dir_cols // 'up.dat')
     call save_matrix(vp(1:nv,1:4), dir_cols // 'vp.dat')
     call save_matrix(hp(1:nh,1:4), dir_cols // 'hp.dat')

end subroutine write_global_cgrid

!***********************************************************************************************************

subroutine gen_global_grid(topo, ph, th, nlon, H_hg, dir_grid, P)!, lonP, latP)

    implicit none
!% receives the topography data
!%
!%  topo(nlon+1, nlon/2+1), lat(nlon/2+1), lon(nlon+1), latp, lonp
!%
!% where topo is a lat x lon grid of height (m) above sea level,
!% lat and lon are the corresponding arbitrarily spaced vectors,
!% and latp and lonp give the position of the rotated pole
!% (all in degrees).
!%
!% Generates whatever grids are necessary, and updates flags as appropriate.
!% Uses the input resolution flags.nph.

     integer, intent(in)::     nlon
     integer			::     nph
     real(dp), allocatable :: ph(:), th(:)
     integer, allocatable   :: topo(:, :)
!     real(dp)     ::     latP, lonP
     type(params) :: P
     real(dp)  ::      dph


     real(dp), allocatable :: ph_ug(:), ta_ug(:), th_ug(:), ph_vg(:), ta_vg(:), th_vg(:)
     real(dp), allocatable :: H_hg(:, :)

     integer ::      statusx, statusy!, statusz
     integer ::      j, js, jn
!     real(dp)  ::      dx, dy, dxR, dyR

     integer ::     nn, ns, nta
     real(dp)  ::      ta_s, ta_n
     integer ::     numLons, numLats

     character(len=*) :: dir_grid

!	real(dp), allocatable :: ph1(:), th1(:), h1(:,:)


!'nph','nta','ta_ug','ta_vg','ph_ug','ph_vg','th_ug','th_vg'
!%======================
!% make some short-cuts
!%======================
nph = P%nph
numLons = nlon+1
numLats = nlon/2+1

    dph=2*pi/nph;    ! the distance between gridpoints

!%==================================================================
!% calculate js (the southernmost index where ocean first appears):
!%==================================================================

js=1;
do while (minval(topo(:, js)) >= 0)
  js=js+1;
end do

!% calculate the corresponding (negative) value of tau:
ta_s=th2ta(th(js), P%coor);

!print *, "th(js), ta_s, th2ta(ta_s)", th(js), ta_s, ta2th(ta_s, P%coor)

!% this will lie in the ns-th coarse tau box south of the equator, where
ns=ceiling(-ta_s/dph);
!%==================================================================
!% calculate jn (the northernmost index where ocean first appears):
!%==================================================================

jn=numLats;
do while (minval(topo(:, jn)) >= 0)
  jn=jn-1;
end do

!% calculate the corresponding value of tau:
ta_n=th2ta(th(jn), P%coor);
!print *, "th, ta_n, th2ta(ta_n)", th(jn), ta_n, th2ta(ta_n, P%coor)

!% this will lie in the nn-th coarse tau box north of the equator, where
nn=ceiling(ta_n/dph);
!print *, "ns", ns
!print *, "nn", nn
!print *, r2d(ta_n), r2d(ta_s)
!%===================================
!% Generate global grid coordinates:
!%===================================
    if (allocated(ph_ug)) deallocate(ph_ug)
    allocate(ph_ug(nph), stat = statusx)
    if (allocated(ph_vg)) deallocate(ph_vg)
    allocate(ph_vg(nph), stat = statusx)

    nta=nn+ns+2

    if (allocated(ta_ug)) deallocate(ta_ug)
    allocate(ta_ug(nta), stat = statusy)
    if (allocated(ta_vg)) deallocate(ta_vg)
    allocate(ta_vg(nta+2), stat = statusy)
    if (allocated(th_ug)) deallocate(th_ug)
    allocate(th_ug(nta), stat = statusy)
    if (allocated(th_vg)) deallocate(th_vg)
    allocate(th_vg(nta+2), stat = statusy)

ph_ug=(/ ( (2*pi*j/nph), j=0,nph-1 ) /)
ph_vg=ph_ug+pi/nph

!ta_vg=(/ ( (j*(2*pi/nph)), j=-ns-2,nn+1 ) /)  ! -ns-1 doesn't always really work...
!th_vg=(/ ( (ta2th(ta_vg(j+ns+2))), j=-ns-2,nn+1 ) /) ! +ns+2 to adjust the lower index
ta_vg=(/ ( ((j-ns-3)*(2*pi/nph)), j=1,nta+2 ) /)  ! -ns-1 doesn't always really work...
!th_vg=(/ ( (ta2th(ta_vg(j), P%coor)), j=1,nn+ns+4 ) /) ! +ns+2 to adjust the lower index
th_vg = ta2th(ta_vg, P%coor)

ta_ug=(/ ( (ta_vg(j)+pi/nph), j=1,nta ) /)
!th_ug=(/ ( (ta2th(ta_ug(j), P%coor)), j=1,nta ) /)
th_ug = ta2th(ta_ug, P%coor)

!%=============================================
!% Interpolate bathymetry to global grid H_hg
!%=============================================
!     call interp_H_grid(numLons, numLats, ph, th, topo, nph, nta, ph_vg, th_ug, H_hg)
	write(*, '("Interpolating topo to ", f5.2, " min resolution global grid: ")', advance = 'no') 360*60/real(nph)
    call bilinear_2d(numLons, numLats, ph, th, real(topo,dp), nph, nta, ph_vg, th_ug, H_hg, P%lib)
	print *, "done."

!%==========================
!% save the coarse C grid
!%==========================
     call save_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call save_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call save_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call save_vector(ta_vg, dir_grid // 'ta_vg.dat')
     call save_vector(th_ug, dir_grid // 'th_ug.dat')
     call save_vector(th_vg, dir_grid // 'th_vg.dat')

!! test bilinear_grid2grid
!	do i = 1, nph
!     H_hg(i,:) = 2*i
!    end do
!	do i = 1, nta
!     H_hg(:,i) = H_hg(:,i)-4*i
!    end do

!	SAVE GENERATED GRID
!     call save_topo(dir_grid // 'topo_orig.dat', nph, nta, ph_vg, th_ug, H_hg)


!    allocate(ph1(nph*10), stat = statusx)
!    allocate(th1(nta*10), stat = statusx)
!    ph1=(/ ( (2*pi*i/(nph*10))+pi/(nph*10), i=0,nph*10-1 ) /)
!    th1=(/ ( ( (th_ug(nta)-th_ug(1))/( nta*10-1 )*i + th_ug(1) ), i=0,nta*10-1 ) /)
!
!      call bilinear_grid2grid(nph, nta, ph_vg, th_ug, real(H_hg,dp), nph*10, nta*10, ph1, th1, h1)!, dir_grid)
!
!      !write to the file
!      call save_topo(dir_grid // 'topo_interp.dat', nph*10, nta*10, ph1, th1, h1)

end subroutine gen_global_grid
!***********************************************************************************************************

subroutine remove_seas(H_hg, nph, nta, Hstar, cmin)

    implicit none
!% receives the topography data H_hg(nph, nta)
!%
!% cmin gives the mimum number of connected points for "sea" to be significant
!% the removed seas are assigned H value Hstar

     integer ::     nta, nph
     real(dp) :: H_hg(nph, nta)
     logical :: in_class(nph, nta)
     integer :: sea_class(nph, nta) ! -1 for land, 0 for unassigned, 1,2,3... for diff seas
     integer :: n_class, cur_class, n_removed
     integer, intent(in) ::      cmin
     real(dp), intent(in) ::      Hstar
     integer, allocatable :: size_class(:)

     integer ::      status
     integer ::      i
     integer :: cph, cta, cph_tmp, cta_tmp

n_class = 0

!ind_b = (H_hg .le. 0 .and. ( (cshift(H_hg, shift=1, dim=1) .ge. 0) .or. &
!                             (cshift(H_hg, shift=-1, dim=1) .ge. 0) .or. &
!                             (cshift(H_hg, shift=1, dim=2) .ge. 0) .or. &
!                             (cshift(H_hg, shift=-1, dim=2) .ge. 0) )
!
!n_b = count(ind_b)

! write(*, '("Removing inland seas/lakes: ")', advance = 'no')

sea_class = 0
cur_class = 0

do cph = 1, nph
     do cta = 1, nta

         if (H_hg(cph,cta) .ge. 0) then
              sea_class(cph,cta) = -1

         elseif (sea_class(cph,cta) .eq. 0) then
!         print *, cph, cta, H_hg(cph,cta), sea_class(cph,cta)
              cur_class = cur_class + 1
              sea_class(cph,cta) = cur_class
              cph_tmp = cph
              cta_tmp = cta
              call mark_connected(H_hg, nph, nta, cph_tmp, cta_tmp, sea_class, cur_class, n_class)
         endif
    enddo
enddo

write(*, '("found ", a, " connected components;")',advance='no') tostring(n_class)

if (allocated(size_class)) deallocate(size_class)
    allocate(size_class(n_class+2), stat = status)

n_removed = -1

do i = -1, n_class

    in_class = sea_class .eq. i
    size_class(i+2) = count(in_class)

    if ((i == -1).or.(size_class(i+2) < cmin)) then
         where (in_class) H_hg = Hstar ! assign height on the land to Hstar
         n_removed = n_removed + 1
    endif

enddo
write(*, '(" removed ", a, " of size < ", a, " cells.")') tostring(n_removed), tostring(cmin)
!write(*,'(200i8)') ( size_class(i+2), i=-1,n_class )

!        print *, " done."

end subroutine remove_seas

!***********************************************************************************************************
subroutine fill_in_const(nph, nta, data, mask)

    implicit none
!% receives mask(nph, nta) and the c2n(nph, nta).
!% mask has the following vals for sponge: -1; water: 0; land: 1.
!% assign data within sponge to the average of the boundary vals

  	 integer, intent(in)		:: mask(nph, nta)
  	 real(wp), intent(inout)	:: data(nph, nta)

     integer ::     nta, nph
     integer :: sponge_class(nph, nta) ! -1 for land, 0 for water, 1,2,3... for diff sponges
     logical :: in_class(nph, nta), water_class(nph, nta),  ind_b(nph, nta)

     integer :: n_class, cur_class, size_class
     integer :: j, status
     integer :: cph, cta, cph_tmp, cta_tmp
     real(wp):: data_avrg_b

n_class = 0
sponge_class = 0
cur_class = 0

do cph = 1, nph
     do cta = 1, nta

         if (mask(cph,cta) == 0) then ! water
              sponge_class(cph,cta) = 0
!         elseif (mask(cph,cta) > 0) then ! land
!              sponge_class(cph,cta) = -1

         elseif (sponge_class(cph,cta) == 0) then ! sponge
!         print *, cph, cta, mask(cph,cta), sponge_class(cph,cta)
              cur_class = cur_class + 1
              sponge_class(cph,cta) = cur_class
              cph_tmp = cph
              cta_tmp = cta
!              call mark_connected(real(mask, kind=wp), nph, nta, cph_tmp, cta_tmp, sponge_class, cur_class, n_class)
			! allow connectivity via land. This should speed things up considerably
              call mark_connected(real(-abs(mask), kind=wp), nph, nta, cph_tmp, cta_tmp, sponge_class, cur_class, n_class)
         endif
    enddo
enddo
! and now remove land from fill-in
where (mask > 0) sponge_class=-1

!write(*, '("found ", a, " connected sponges;")',advance='no') tostring(n_class)

water_class = sponge_class .eq. 0
do j = 1, n_class

    in_class = (sponge_class .eq. j)
    size_class = count(in_class)

    if (size_class>nph/10) then ! fill in only significantly large sponge regions; smaller ones are efficiently treated later
		ind_b = (water_class .and. ( cshift(in_class, shift=1, dim=1) .or. &
		                             cshift(in_class, shift=-1, dim=1).or. &
		                             cshift(in_class, shift=1, dim=2).or. &
		                             cshift(in_class, shift=-1, dim=2) ))

		data_avrg_b = sum(data, mask=ind_b)/count(ind_b) ! calc the average of data over the sponge-water boundary
		where (in_class) data = data_avrg_b ! assign data within sponge to the average of the boundary vals
	endif
enddo

end subroutine fill_in_const

!***********************************************************************************************************
subroutine mark_connected(H_hg, nph, nta, cph, cta, sea_class, cur_class, n_class)
    implicit none

     integer ::     nta, nph
     real(dp) :: H_hg(nph, nta)
     integer :: sea_class(nph, nta) ! -1 for land, 0 for unassigned, 1,2,3... for diff seas
     integer :: cur_class, n_class
     integer :: cph, cta, cphl, cphr
     integer :: buf(nph*nta*3/4, 2) ! FIFO buffer for neighbours' coordinates
     integer :: bufp, bufe

n_class = n_class + 1

!write( *, '(4i6)', advance ='no' )  cph, cta, H_hg(cph,cta), n_class,cur_class

bufp = 1
bufe = 1

call push(buf, cph, cta, bufe)

do while (pop(buf, cph, cta, bufe, bufp))

     if (cta < nta) then
       if ((H_hg(cph, cta+1) < 0) .and. (sea_class(cph, cta+1) .eq. 0)) then
         !write( *, '(" 1")', advance = 'no')
         sea_class(cph, cta+1) = cur_class
         !print *, " 1", cph, cta+1
         call push(buf, cph, cta+1, bufe)
       endif
     endif

     if (cta > 1) then
       if ((H_hg(cph,cta-1) < 0) .and. (sea_class(cph,cta-1) .eq. 0)) then
!         write( *, '(" 4")', advance = 'no')
         sea_class(cph,cta-1) = cur_class
!         print *, " 4", cph, cta-1
         call push(buf, cph, cta-1, bufe)
       endif
     endif

     cphl=1 + modulo((cph-1)-1,nph)
     cphr=1 + modulo((cph+1)-1,nph)
!     print *, cphl, cphr

     if ((H_hg(cphr,cta) < 0) .and. (sea_class(cphr,cta) .eq. 0)) then
 !    print *, " 2"
         sea_class(cphr,cta) = cur_class
         call push(buf, cphr, cta, bufe)
     endif

     if ((H_hg(cphl,cta) < 0) .and. (sea_class(cphl,cta) .eq. 0)) then
     !print *, " 3"
         sea_class(cphl,cta) = cur_class
         call push(buf, cphl, cta, bufe)
     endif

enddo

!print *, bufe, bufp

end subroutine mark_connected

!***********************************************************************************************************

subroutine push(buf, cph, cta, bufe)
    implicit none

     integer :: cph, cta
     integer :: buf(:,:) !(nph*nta*2/3, 2)
     integer :: bufe

buf(bufe, 1) = cph
buf(bufe, 2) = cta
bufe = bufe + 1

end subroutine push

logical function pop(buf, cph, cta, bufe, bufp)
    implicit none

     integer :: cph, cta
     integer :: buf(:,:) !(nph*nta*2/3, 2)
     integer :: bufe, bufp
if (bufe == bufp) then
     pop = .false.
     return
else
     cph = buf(bufp, 1)
     cta = buf(bufp, 2)
     bufp = bufp + 1
     pop = .true.
     return
endif

end function pop

!***********************************************************************************************************

subroutine remove_bays(H_hg, nph, nta, Hstar)!,dir_grid)

    implicit none
!% receives the topography data H_hg(nph, nta) and removes bays
!% and inlets which are just one grid cell wide.
!% the removed seas are assigned H value Hstar

     integer, intent(in)  ::     nta, nph
     real(dp)			  :: H_hg(nph, nta)
     real(dp), intent(in) :: Hstar
     ! true if boundary is present on the  left, right, above, below
     logical :: l_b(nph, nta), r_b(nph, nta), u_b(nph, nta), b_b(nph, nta)
!     integer :: tmp(nph, nta)
     logical :: ind_3b(nph, nta) ! true when 3 boundaries are around a water cell
     integer :: n_removed

!	character(len=*) :: dir_grid

write(*, '("Removing bays and inlets: ")', advance='no')

n_removed = 1
do while (n_removed > 0)

! find position of cells which has 3 sides on the boundary
r_b = (H_hg < 0.) .and. (cshift(H_hg, shift=1, dim=1) >= 0.)
l_b = (H_hg < 0.) .and. (cshift(H_hg, shift=-1, dim=1) >= 0.)
b_b = (H_hg < 0.) .and. (cshift(H_hg, shift=1, dim=2) >= 0.)
u_b = (H_hg < 0.) .and. (cshift(H_hg, shift=-1, dim=2) >= 0.)

ind_3b = (l_b .and. r_b .and. u_b) .or. (l_b .and. r_b .and. b_b) .or. &
         (l_b .and. b_b .and. u_b) .or. (r_b .and. b_b .and. u_b)

    where (ind_3b) H_hg = Hstar ! assign height on the land to Hstar

    n_removed = count(ind_3b)
    write(*, '(a)', advance = 'no') tostring(n_removed)
    if (n_removed == 0) then
    	write(*, '(".")')
    else
    	write(*, '(" + ")', advance = 'no')
    endif

end do

!call save_matrix(cshift(H_hg, shift=1, dim=1), dir_grid // 'l_shift.dat')
!
!tmp = 0
!where (l_b) tmp=1
!call save_matrix(tmp, dir_grid // 'l_b.dat')
!tmp = 0
!where (r_b) tmp=1
!call save_matrix(tmp, dir_grid // 'r_b.dat')
!tmp = 0
!where (u_b) tmp=1
!call save_matrix(tmp, dir_grid // 'u_b.dat')
!tmp = 0
!where (b_b) tmp=1
!call save_matrix(tmp, dir_grid // 'b_b.dat')

!        print *, "Done."

end subroutine remove_bays

!***********************************************************************************************************

subroutine adjust_min_depth(H_hg, nph, nta, P, dir_grid)

    implicit none
!% receives the topography data H_hg(nph, nta) and
!% adjusts the H_hg by setting the min depth
! scheme 1: a simple minimum depth
! scheme 2: min number of Grid Points Per Wavelength

     integer               ::     nta, nph
     real(dp)              :: H_hg(nph, nta)
     real(dp)              :: ta_ug(nta)

     type(params), intent(in)	:: P
     type (tide_params)	   		:: t_params(P%ncpts)

     integer   :: scheme, gppw, coor
     real(dp)  ::      hmin
     real(dp)  ::      g                 ! surface gravity, m/s^2
     real(dp)  ::      re                ! average earth's radius, m
!     real(dp)  ::      omega             ! angular velocity, rad/s

     real(dp) :: dx(nph, nta), kmax(nph, nta), hmin_m(nph, nta)
     integer  :: i
     real(dp) :: omegamax

     character(len=*) :: dir_grid

coor = P%coor
scheme = P%hmin_scheme
hmin =  P%hmin
gppw = P%gppw
g = P%g
re = P%re
!omega = P%omega

if (scheme == 1) then

     where ((H_hg < 0.).and.(H_hg > -hmin)) H_hg = -hmin

elseif (scheme == 2) then

!    omegamax = 2 * omega ! assume no higher constituents (m2 has the highest frequency)
	t_params = get_params(P%ncpts, P%cpts)
	omegamax = maxval(t_params%omega0)

     if (coor == 1) then ! Mercator

          call load_vector(ta_ug, dir_grid // 'ta_ug.dat')

          do i = 1, nph
              dx(i, :) = (2*pi*re/nph)*(1/cosh(ta_ug))
          enddo

     else if  (coor == 2) then ! lat-lon

          dx = 2*pi*re/nph

     endif

  kmax = 2*pi/(gppw*dx)   ! from min_wavelength = 2*pi/kmax = *MUST BE* = dx*gppw;

  hmin_m = g * (omegamax/kmax)**2

  where (hmin_m < hmin) hmin_m = hmin ! adjust for global minimum

  where ((H_hg > -hmin_m).and.(H_hg < 0.)) H_hg = -hmin_m ! now apply to H

endif

end subroutine adjust_min_depth
!**********************************************************************************************************

subroutine calc_sponge_prefactor(hfactor, omega0, P, nh, dir_cols, dir_grid)

implicit none

	type(params) :: P
    integer, intent(in)    :: nh
	real(wp)           :: omega0

     real(wp), allocatable	:: ph_vg(:), th_ug(:), ph_h(:), th_h(:)
     real(wp), allocatable	:: ph_h_nonrot(:), th_h_nonrot(:)
	 real(wp)           :: th0, ph0
	 integer, allocatable   :: hp(:,:)

     integer            :: j, istat
     character(len = *) :: dir_grid, dir_cols!, dir_mats, N_data_dir

     logical, allocatable, dimension(:) :: hfactor


allocate(hfactor(nh), stat = istat)
hfactor = 0

! load hp, ph_vg, ta_h
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')

    allocate(ph_h(nh), th_h(nh), stat = istat)
! and fill them
	ph_h = (/ ( (ph_vg(hp(j, 1))), j=1,nh ) /)
	th_h = (/ ( (th_ug(hp(j, 2))), j=1,nh ) /)
	deallocate(ph_vg, th_ug)

! poles
th0 = d2r(P%latP)
ph0 = d2r(P%lonP)
! transform to the non-rotated coords
	call calc_latlon_nonrot(ph_h, th_h, ph0, th0, ph_h_nonrot, th_h_nonrot)

	where ((abs(th_h_nonrot)) < asin(.5*omega0/P%omega))  hfactor = .true.

	deallocate(ph_h, th_h, ph_h_nonrot, th_h_nonrot)

end subroutine calc_sponge_prefactor
!***********************************************************************************************************

subroutine generate_sponge(cpts, P, GD, dir_grid, dir_cols, dir_mats)

    implicit none
!% Receives cn on the h-grid. Phase speed cn is calculated with:
							! cn_scheme: 1 -- WKB; 2 -- evaluate Chebyshev series for cn(H)
!% Marks the sponge layer where omega/cn < kmax (kmax -- maximum resolved wavenumber)
							! min number of Grid Points Per Wavelength is given by P%igw_gppw
     type(params), intent(in)	:: P
     type(grid_dims), intent(in):: GD
     character(len=*)			:: cpts

     integer				:: nu, nv, nh, nta, nph
     character	:: lib
     real(dp), allocatable	:: ph_vg(:), th_ug(:), ta_ug(:), H_hg(:, :), H_hg_orig(:,:)
     real(dp), allocatable	:: c2n_interp(:,:),c2n_orig(:,:), c2n_hg(:,:), c2n_hg_orig(:,:)
     real(dp), allocatable	:: c2n_u(:,:), c2n_v(:,:), c2n_h(:,:), c2n_u_tmp(:,:), c2n_v_tmp(:,:)!, dc2n_h(:,:)

     integer, allocatable 	:: hp(:, :), up(:,:), vp(:,:)
     real(dp), allocatable	:: H_h(:), check_insponge_seas(:,:) !, H_u(:), H_v(:)

     logical, allocatable, dimension(:) :: hfactor
     integer, allocatable	:: issponge_hg(:, :, :)
	 real(wp), allocatable	:: issponge_u(:, :), issponge_v(:, :), issponge_h(:, :), issponge_hg_sm(:, :, :)
	 integer				:: edge, cedge
!     integer				:: lb, rb, bb, ub
  	 logical, allocatable	:: ind_b(:,:)
     real(dp), allocatable	:: ph_b(:), th_b(:)
!	 real(dp)				:: sm_dist
	integer					:: switch ! switch on/off proper calcs on neighbours for bilinear interpolation

	real(dp), allocatable	:: N(:), z(:)
	real(dp)				:: Nbar

     type (tide_params)	   	:: t_params
     integer				:: tide_types(P%ncpts), ntypes, ctypes

     integer   :: gppw, coor, cn_scheme, nmodes, nz, N_hor_var, n_b, c_b
     real(dp)  ::      re                ! average earth's radius, m
!     real(dp)  ::      omega             ! angular velocity, rad/s

     real(dp) :: dx(GD%nph, GD%nta), kmax(GD%nph, GD%nta), hmin_m(GD%nph, GD%nta)
     integer  :: j, k, l, istat, loc, dims(2)

     real(dp) :: omegamax

!     logical :: l_b(GD%nph, GD%nta), r_b(GD%nph, GD%nta), u_b(GD%nph, GD%nta), b_b(GD%nph, GD%nta)
!     logical :: ind_3b(GD%nph, GD%nta) ! true when 3 boundaries are around a water cell
	 type (triplet) 		:: h2v, h2u


     character(len=*) :: dir_grid, dir_cols, dir_mats
!------------------------------------------------------------------------------------------------
!	SHORTCUTS
coor = P%coor
nmodes = P%n_modes
cn_scheme = P%cn_scheme
N_hor_var=P%N_hor_var
gppw = P%igw_gppw
re = P%re
lib = P%lib

nph = GD%nph
nta = GD%nta
nu = GD%nu
nv = GD%nv
nh = GD%nh
!------------------------------------------------------------------------------------------------
write(*, '("Generating sponge layer(s) for ",i2," modes. ")') nmodes

	if (P%smooth_type==0) then
		call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
		call load_alloc_vector(H_h, dir_cols // 'H_h.dat')
	else
		call load_alloc_matrix(H_hg, dir_grid // 'H_sht_hg.dat')
		call load_alloc_vector(H_h, dir_cols // 'H_sht_h.dat')
		call load_alloc_matrix(H_hg_orig, dir_grid // 'H_hg.dat')
	endif

	call load_alloc_matrix(hp, dir_cols // 'hp.dat')
	call load_alloc_matrix(up, dir_cols // 'up.dat')
	call load_alloc_matrix(vp, dir_cols // 'vp.dat')

	allocate(issponge_hg(nph, nta, nmodes),issponge_hg_sm(nph, nta, nmodes), &
			 issponge_u(nu, nmodes), issponge_v(nv, nmodes), issponge_h(nh, nmodes), stat=istat)
	allocate(c2n_hg(nph, nta), stat=istat)
	if ( (N_hor_var==1).and.(P%smooth_type >= 1) ) then
		allocate(c2n_hg_orig(nph, nta), stat=istat)
	endif

!	Generate TWO different sponges for semidiurnal and diurnal constituents
	! check which are present
	tide_types = 0
!	write(*,'(a)') cpts
	do j = 1, P%ncpts
		if ( cpts(2*j:2*j) .eq. '1') then ! diurnal
			tide_types(j) = 1
		elseif ( cpts(2*j:2*j) .eq. '2') then ! semidiurnal
			tide_types(j) = 2
		endif
	enddo
	call unique(tide_types, ntypes)

     if (coor == 1) then ! Mercator
          call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
          do j = 1, nph
              dx(j, :) = (2*pi*re/nph)*(1/cosh(ta_ug))
          enddo
     else if  (coor == 2) then ! lat-lon
          dx = 2*pi*re/nph
     endif

  kmax = 2*pi/(gppw*dx)   ! from min_wavelength = 2*pi/kmax = *MUST BE* = dx*gppw;
!------------------------------------------------------------------------------------------------

do ctypes = 1, ntypes ! a tiny cycle to make separate sponge for diurnal and semidiurnal tides

	issponge_hg = 0
	issponge_hg_sm = 0
	issponge_u = 0
	issponge_v = 0
	issponge_h = 0

	write(*, '("")')
	if (tide_types(ctypes) == 1) then
		t_params = get_pars('k1') ! for all diurnal use 'k1'
		write(*,'(a, e10.4,a)') 'Diurnal frequencies (',t_params%omega0, ')'
	elseif (tide_types(ctypes) == 2) then
		t_params = get_pars('m2') ! for all semidiurnal use 'm2'
		write(*,'(a, e10.4, a)') 'Semidiurnal frequencies (', t_params%omega0, ')'
	endif

!	t_params = get_params(ncpts, cpts)
!	omegamax = maxval(t_params%omega0)
	omegamax = t_params%omega0 ! name omagemax for historical reasons

	!------------------------------------------------------------------------------------------------
	if (cn_scheme == 1) then ! WKB formula for cn(H)

		! Calculate N_bar (i.e. z average of N)
		call load_alloc_vector(z, dir_grid // 'z.dat')
		nz = size(z)
		call load_alloc_vector(N, dir_grid // 'N_avrg.dat')
		Nbar = 1/(z(nz)-z(1)) * sum( (N(2:nz)+N(1:nz-1))/2 * (z(2:nz)-z(1:nz-1)) ) ! trapz rule averaging

		if (allocated(c2n_interp)) deallocate(c2n_interp)
		allocate(c2n_interp(nh, nmodes), stat=istat)
	  	do k = 1, nmodes
	  		c2n_interp(:, k) = (Nbar*H_h/pi/k)**2
			hmin_m = k * pi/Nbar * (omegamax/kmax)
	!		call save_matrix(hmin_m, dir_mats//'temp/'//'hmin_'//tostring(k)//'.dat')

		    if (P%smooth_type==0) then
		    	where ((H_hg > 0.).and.(H_hg < hmin_m)) issponge_hg(:,:,k) = 1
		    else
		    	where ( ((H_hg > 0.).and.(H_hg < hmin_m)) .or. &
		    		  ((H_hg_orig > 0.).and.(H_hg_orig < hmin_m)) ) issponge_hg(:,:,k) = 1
		    endif

			write(*, '("Minimal resolved depths for cpt=*", i1, ": ")', advance='no'), tide_types(ctypes)
		    write(*, '("mode ",i2,": ",f6.1,"m; ")', advance='no') k, minval(hmin_m(:, nta/2))
	!	    print *, hmin_m(1,1)
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!***************	IMPORTANT: ADD WKB, VARIABLE N SPONGE GENERATION ******************************
		! ps and correct case (N_hor_var==0) to be like (N_hor_var==1)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	elseif (cn_scheme == 2) then !

		! load cn
		call load_alloc_matrix(c2n_interp, dir_grid // 'itm/' // 'c2n_h.dat')

		if (N_hor_var==0) then ! use cn(H) to find hmin
			do k = 1, nmodes
				! find H_h at which the value of c2n is the closest to critical
			     if (coor == 1) then ! Mercator
			          do j = 1, nta ! kmax depends on lat only
			          	loc = minloc(c2n_interp(:, k) - omegamax/kmax(1, j), 1, mask = c2n_interp(:, k) > omegamax/kmax(1, j))
			          	hmin_m(:, j) = H_h(loc)
			          enddo
			     else if  (coor == 2) then ! lat-lon, dx is constant
			          loc = minloc(c2n_interp(:, k) - omegamax/kmax(1, 1), 1, mask = c2n_interp(:, k) > omegamax/kmax(1, 1))
			          if (loc == 0) then
			          	write(*,*) ""
			            write(*, '("ERROR: Mode ",i2," cannot be resolved accurately within the given criteria")') k
			          	write(*,*) 'Consider increasing resolution or reducing the number of resolved IT modes'
			          	stop
			          endif
			          hmin_m = H_h(loc)
		!	          print *, H_h(loc)
			     endif

			    if (P%smooth_type==0) then
			    	where ((H_hg > 0.).and.(H_hg < hmin_m)) issponge_hg(:,:,k) = 1
			    else
			    	where ( ((H_hg > 0.).and.(H_hg < hmin_m)) .or. &
			    		  ((H_hg_orig > 0.).and.(H_hg_orig < hmin_m)) ) issponge_hg(:,:,k) = 1
			    endif

			    if (k==1) then
			    	write(*, '("Minimal resolved depths for cpt=*", i1, ": ")', advance='no'), tide_types(ctypes)
			    endif
			    write(*, '("mode ",i2,": ",f6.1,"m; ")', advance='no') k, minval(hmin_m(:, nta/2))
		!	    print *, '(Before removing sponge isles) Mode', k, 'Relative sponge area %', 100*count(issponge_hg(:,:,k)==1)/count(H_hg>0)
			enddo

		elseif (N_hor_var==1) then ! use cn(H,N) to find hmin

			if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
				call load_alloc_matrix(c2n_orig, dir_grid // 'itm/' // 'c2n_h_orig.dat')
			endif

			do k = 1, nmodes

				! convert cn_h(:,k) into a sparse matrix
				c2n_hg = 0
				do j = 1, nh
					c2n_hg(hp(j,1),hp(j,2)) = c2n_interp(j, k)
				end do
				if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
					c2n_hg_orig = 0
					do j = 1, nh
						c2n_hg_orig(hp(j,1),hp(j,2)) = c2n_orig(j, k)
					end do
				endif


			    if ((P%smooth_type==0).or.(P%baro_on_smoothed==1)) then
			    	where ((H_hg > 0.).and.(c2n_hg <= omegamax/kmax)) issponge_hg(:,:,k) = 1
			    else
			    	where ( ((H_hg > 0.).and.(c2n_hg <= omegamax/kmax)) .or. &
			    		    ((H_hg_orig > 0.).and.(c2n_hg_orig <= omegamax/kmax)) ) issponge_hg(:,:,k) = 1
			    endif

			       if ( count(issponge_hg(:,:,k)==1)/count(H_hg>0) > .9 ) then
			         	write(*,*) ""
			            write(*, '("ERROR: Mode ",i2," cannot be resolved accurately within the given criteria")') k
			          	write(*,*) 'Consider increasing resolution or reducing the number of resolved IT modes'
		!	    print *, '(Before removing sponge isles) Mode', k, 'Relative sponge area %', 100*count(issponge_hg(:,:,k)==1)/count(H_hg>0)
			          	stop
			       endif
	!		    write(*, '("mode ",i2,": ",f6.1,"m; ")', advance='no') k, minval(H_hg, ((H_hg > 0.).and.(issponge_hg(:,:,k) == 0)) )
			enddo
		else
			write(*,*) "ERROR in subroutine generate_sponge: parameter N_hor_var must be 0 or 1"
			stop
		endif

	endif
		write(*,*) ""

	!%==============================
	!% remove water-islands in the sponge
	!%==============================
		allocate(check_insponge_seas(nph, nta), stat=istat)
		do k = 1, nmodes
			check_insponge_seas = issponge_hg(:,:,k) - .5 ! so that sponge points have vals .5 and off-sponge have vals -.5
			where (H_hg <= 0.) check_insponge_seas = 1._dp ! so that true land is not considered


			write(*, '("Removing in-sponge seas/lakes: ")', advance = 'no')
			call remove_seas(check_insponge_seas, nph, nta, .5_dp, 2*nph)  ! new depth is 0.5, 2*nph is alot (to remove large inland seas!)
			where ((H_hg > 0.).and.(check_insponge_seas > 0)) issponge_hg(:,:,k) = 1
		enddo
		deallocate(check_insponge_seas) ! will be needed

		write(*, '("Relative sponge area for cpt=*", i1, ": ")', advance='no'), tide_types(ctypes)
		do k = 1, nmodes
			write(*, '("mode ",i2,": ", f4.1,"%; ")', advance='no') k, real(100*count(issponge_hg(:,:,k)==1))/count(H_hg>0)
		enddo

		call save_matrix_3d(issponge_hg, dir_grid// 'itm/' //'issponge_hg'//tostring(tide_types(ctypes))//'.dat')

		write(*,*) ""


	!%========================
	!% Adjust cn in the sponge
	!%========================
	write(*, '("Smooth cn within the sponge area for cpt=*", i1, ": ")', advance='no'), tide_types(ctypes)
	call load_alloc_sparse(h2u, dir_mats // 'h2u.dat')
	call load_alloc_sparse(h2v, dir_mats // 'h2v.dat')

	call load_alloc_matrix(c2n_u, dir_grid // 'itm/' // 'c2n_u.dat')
	call load_alloc_matrix(c2n_v, dir_grid // 'itm/' // 'c2n_v.dat')
	dims = shape(c2n_u)
	allocate(c2n_u_tmp(dims(1),dims(2)), stat=istat)
	dims = shape(c2n_v)
	allocate(c2n_v_tmp(dims(1),dims(2)), stat=istat)
!	call load_alloc_matrix(c2n_h, dir_grid // 'itm/' // 'c2n_h.dat')
!	call load_alloc_matrix(dc2n_h, dir_grid // 'itm/' // 'dc2n_h.dat')

	do k = 1, nmodes
		write(*, '("mode ",i2,"; ")', advance='no') k
		! convert cn_h(:,k) into a sparse matrix
		c2n_hg = 0.
		do j = 1, nh
			c2n_hg(hp(j,1),hp(j,2)) = c2n_interp(j, k)
		end do

		if (P%sponge_smooth==1) then
			! use front-fill to initalize c2n inside the sponge (faster, cruder)
			call front_fill_in(c2n_hg, issponge_hg(:,:,k))
		elseif (P%sponge_smooth==2) then
			! before smoothing initialize c2n in the sponge with a constant (average of boundary vals)
			! sponge: -1; water: 0; land: 1
			call fill_in_const(nph, nta, c2n_hg, log2int((H_hg <= 0.)) - issponge_hg(:,:,k))
			! relaxation smoothing (solve laplace eq within sponge)
			call laplace_relax(c2n_hg, issponge_hg(:,:,k))
		endif

		! now save adjusted c2n
		! h-grid
		do j = 1, nh
			c2n_interp(j, k) = c2n_hg(hp(j,1),hp(j,2))
		enddo

		! u-grid
		call coo_vec_mul(h2u, c2n_interp(:, k), P%lib, c2n_u_tmp(:, k))
		do j = 1, nu
			if (issponge_hg(up(j,1),up(j,2),k)==1) then
				c2n_u(j, k) = c2n_u_tmp(j, k)
			endif
		end do

		! v-grid
		call coo_vec_mul(h2v, c2n_interp(:, k), P%lib, c2n_v_tmp(:, k))
		do j = 1, nv
			if (issponge_hg(vp(j,1),vp(j,2),k)==1) then
				c2n_v(j, k) = c2n_v_tmp(j, k)
			endif
		end do

	enddo
	if (P%sponge_smooth==0) then
		write(*, '("none!")')
	else
		write(*, '("done.")')
	endif

	call save_matrix(c2n_interp, dir_grid // 'itm/' // 'c2n_h'//tostring(tide_types(ctypes))//'.dat')
	call save_matrix(c2n_u, dir_grid // 'itm/' // 'c2n_u'//tostring(tide_types(ctypes))//'.dat')
	call save_matrix(c2n_v, dir_grid // 'itm/' // 'c2n_v'//tostring(tide_types(ctypes))//'.dat')

	deallocate(c2n_u, c2n_v, c2n_u_tmp, c2n_v_tmp)

	!%============================================================
	!% smooth sponge to make it more gradual
	!%============================================================
	edge=3*gppw
! SWITCH BETWEEN CALCULATING DISTANCES TO THE SPONGE-OCEAN BOUNDARY DIRECTLY (VERY SLOW, NOT RECOMMENDED)
!			 AND A LOOP TO CALCULATE GRID POINTS NOT MORE THAN edge AWAY FROM THE SPONGE-OCEAN BOUNDARY
switch = 0 ! fast version, calc # of gridcells inbetween
!switch = 1 ! calculate actual geodesic distances and find the minimum

if (switch == 0) then

	allocate (ind_b(nph, nta), stat=istat)

	do k = 1, nmodes

		c2n_hg = 0.
		do j = 1, nh
			c2n_hg(hp(j,1),hp(j,2)) = c2n_interp(j, k)
		end do

	! Mark issponge_hg_sm based on how many shifting iterations it is from the ocean

	! zero-level iteration is the 'ocean' boundary of the sponge (where issponge==0 and cn/=0)
		ind_b = ( (c2n_hg > 0) .and. (issponge_hg(:,:,k) == 0) .and. &
		                           ( (cshift(issponge_hg(:,:,k), shift=1, dim=1) == 1) .or. &
		                             (cshift(issponge_hg(:,:,k), shift=-1, dim=1) == 1) .or. &
		                             (cshift(issponge_hg(:,:,k), shift=1, dim=2) == 1) .or. &
		                             (cshift(issponge_hg(:,:,k), shift=-1, dim=2) == 1) ) )
		! loop up to #edge shifts away from the boundary
		do cedge = 1,edge-1
			where ((issponge_hg(:,:,k) == 1).and.(.not.(ind_b)).and.&
		                           ( (cshift(ind_b, shift=1, dim=1)) .or. &
		                             (cshift(ind_b, shift=-1, dim=1)) .or. &
		                             (cshift(ind_b, shift=1, dim=2)) .or. &
		                             (cshift(ind_b, shift=-1, dim=2)) ) )
				issponge_hg_sm(:,:,k) = real(cedge, kind=wp)/edge
				ind_b = .true.
			end where
		enddo

		! fill in value for the bulk of the sponge away from the boundary
		where ((issponge_hg(:,:,k) == 1) .and. (issponge_hg_sm(:,:,k)==0)) issponge_hg_sm(:,:,k) = 1
	enddo

	deallocate(ind_b)
elseif (switch == 1) then
!	A BIT WHERE I CALCULATE DISTANCE TO THE SPONGE-OCEAN BOUNDARY DIRECTLY (VERY SLOW, DITCHED IT)
	call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
	call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
	dx = dx/re ! adjust dx to represent distances on a unit sphere

	allocate (ind_b(nph, nta), stat=istat)

	do k = 1, nmodes

		c2n_hg = 0.
		do j = 1, nh
			c2n_hg(hp(j,1),hp(j,2)) = c2n_interp(j, k)
		end do

	! the 'ocean' boundary of the sponge (where issponge==0 and cn/=0)
		ind_b = ( (c2n_hg > 0) .and. (issponge_hg(:,:,k) == 0) .and. &
		                           ( (cshift(issponge_hg(:,:,k), shift=1, dim=1) == 1) .or. &
		                             (cshift(issponge_hg(:,:,k), shift=-1, dim=1) == 1) .or. &
		                             (cshift(issponge_hg(:,:,k), shift=1, dim=2) == 1) .or. &
		                             (cshift(issponge_hg(:,:,k), shift=-1, dim=2) == 1) ) )
!		! ocean boundary within the 'close-by' box
		n_b = count(ind_b)
	! lat-lon coords of the boundary
		if (allocated(ph_b)) deallocate(ph_b)
		if (allocated(th_b)) deallocate(th_b)
		allocate (ph_b(n_b),th_b(n_b), stat=istat)

		c_b=0
		do j = 1, nph
		    do l = 1, nta
				if (ind_b(j,l)) then
					c_b = c_b + 1
					ph_b(c_b) = ph_vg(j)
					th_b(c_b) = th_ug(l)
				endif
			enddo
		enddo
		! for every sponge point calc distance to the ocean boundary
		do j = 1, nph
		    do l = 1, nta
				if (issponge_hg(j,l,k) == 1) then
					issponge_hg_sm(j,l,k) = minval( dist_sphere(ph_vg(j:j), th_ug(l:l), ph_b, th_b) ) / (edge*dx(j,l))
				endif
			enddo
		enddo
	enddo
	! where distance to the boundary is greater than edge*dx(j,l) level off issponge_hg_sm to 1
	where (issponge_hg_sm > 1) issponge_hg_sm = 1

	deallocate(ind_b, ph_b, th_b)
endif
	deallocate(c2n_interp)

	! TEMPRORARY, INVESTIGATE!
	! Smooth sponge regions unphysically affect coastally trapped waves.
	! For costally trapped waves (K1) this may result in regions of negative conversion (hence, bad convergence)
	! Increasing sponge_damp by a factor of 10 removes those. Hence, an ad-hoc fix has been implemented, that multiplies the damping coefficient by 10 for CTW
	if (P%ctw_fix==1) then
		call calc_sponge_prefactor(hfactor, omegamax, P, nh, dir_cols, dir_grid)
		do j = 1, nh
			if (.not. hfactor(j)) then
				issponge_hg_sm(hp(j,1),hp(j,2),:) = issponge_hg(hp(j,1),hp(j,2),:) ! Introduce the abrupt sponge past the critical latitude.
!				issponge_hg_sm(hp(j,1),hp(j,2),:) = 10*issponge_hg_sm(hp(j,1),hp(j,2),:)
			endif
		end do
	endif

!	call save_matrix_3d(issponge_hg_sm, dir_grid// 'itm/' //'issponge_hg_sm'//tostring(tide_types(ctypes))//'.dat')

	!	Save issponge for column-vectors
	do j = 1, nmodes
	!	now map back into u,v,h cols
	    do k = 1, nh
	!		issponge_h(k, j) = issponge_hg(hp(k,1),hp(k,2), j)
			issponge_h(k, j) = issponge_hg_sm(hp(k,1),hp(k,2), j)
		end do
	    do k = 1, nu
	!		issponge_u(k, j) = issponge_hg(up(k,1),up(k,2), j)
			issponge_u(k, j) = issponge_hg_sm(up(k,1),up(k,2), j)
		end do
	    do k = 1, nv
	!		issponge_v(k, j) = issponge_hg(vp(k,1),vp(k,2), j)
			issponge_v(k, j) = issponge_hg_sm(vp(k,1),vp(k,2), j)
		end do
	end do

		call save_matrix(issponge_u, dir_grid//'itm/'//'issponge_u'//tostring(tide_types(ctypes))//'.dat')
		call save_matrix(issponge_v, dir_grid//'itm/'//'issponge_v'//tostring(tide_types(ctypes))//'.dat')
		call save_matrix(issponge_h, dir_grid//'itm/'//'issponge_h'//tostring(tide_types(ctypes))//'.dat')

enddo ! end the tiny cycle to make separate sponge for diurnal and semidiurnal tides

deallocate(H_h, H_hg, c2n_hg)
if ((P%smooth_type>=1).and.(P%baro_on_smoothed/=1)) then
	deallocate(H_hg_orig)
	if (N_hor_var==1) then
		deallocate(c2n_orig, c2n_hg_orig)
	endif
endif

deallocate(issponge_hg, issponge_hg_sm, issponge_u, issponge_v, issponge_h)

end subroutine generate_sponge
!***********************************************************************************************************

!********************************************************************
! Write the parameters of the grid into a specified file
!********************************************************************
subroutine write_GD(file_name, GD)

  implicit none

  character(len=*) :: file_name
  type(grid_dims):: GD

  integer :: fh = 15

  open(fh, file=file_name, action='write', status='replace')

  write(fh, '(a)') '! Basic parameters of the generated grid'
  write(fh, '(a)') 'coor = ' // tostring(GD%coor)
  write(fh, '(a)') 'nph = ' // tostring(GD%nph)
  write(fh, '(a)') 'nta = ' // tostring(GD%nta)
  write(fh, '(a)') 'nu = ' // tostring(GD%nu)
  write(fh, '(a)') 'nv = ' // tostring(GD%nv)
  write(fh, '(a)') 'nh = ' // tostring(GD%nh)
  write(fh, '(a)') 'np = ' // tostring(GD%np)
  write(fh, '(a)') 'lonp = ' // tostring(GD%lonP)
  write(fh, '(a)') 'latp = ' // tostring(GD%latP)

  close(fh)

end subroutine

!********************************************************************
! Compare if two GD structures are identical
!********************************************************************
logical function cmpr_GD(GD, GD2)

  implicit none

  type(grid_dims), intent(in) :: GD, GD2
  logical :: cmpr_GD

 cmpr_GD = 	(GD%coor == GD2%coor ).and.(GD%nph == GD2%nph ).and. &
 			(GD%nta == GD2%nta ).and.(GD%nu == GD2%nu ).and. &
 			(GD%nv == GD2%nv ).and.(GD%nh == GD2%nh ).and. &
 			(GD%np == GD2%np ).and.(GD%lonP == GD2%lonP ).and.(GD%latP == GD2%latP)

end function


!********************************************************************
! Read the parameters of the grid into a specified file
!********************************************************************
subroutine read_GD(file_name, GD)

  implicit none

  character(len=*) :: file_name
  type(grid_dims):: GD

  character(len=100) :: buffer, label
  integer :: pos, pos_
  integer :: fh
  integer :: ios
  integer :: line

fh = 15; ios = 0; line = 0;

  open(fh, file=file_name, action='read')

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos_ = scan(buffer, '!')
        pos = scan(buffer, '=')
        !print *, pos_, pos
        if ( (pos>0).and.((pos<pos_).or.(pos_==0)) ) then

	        label = buffer(1:pos-1)
	        if (pos_ > 0) then
	        	buffer = trim(buffer(pos+1:pos_-1))
	        else
	        	buffer = trim(buffer(pos+1:))
	        end if

	        select case (label)
!********************************************************************
	        case ('coor')
	           read(buffer, *, iostat=ios) GD%coor
	           !print*, 'Read coor: ', GD%coor
	        case ('nph')
	           read(buffer, *, iostat=ios) GD%nph
	           !print*, 'Read nph: ', GD%nph
	        case ('nta')
	           read(buffer, *, iostat=ios) GD%nta
	           !print*, 'Read nta: ', GD%nta
	        case ('nu')
	           read(buffer, *, iostat=ios) GD%nu
	           !print*, 'Read nu: ', GD%nu
	        case ('nv')
	           read(buffer, *, iostat=ios) GD%nv
	           !print*, 'Read nv: ', GD%nv
	        case ('nh')
	           read(buffer, *, iostat=ios) GD%nh
	           !print*, 'Read nh: ', GD%nh
	        case ('np')
	           read(buffer, *, iostat=ios) GD%np
	           !print*, 'Read np: ', GD%np
	        case ('lonp')
	           read(buffer, *, iostat=ios) GD%lonP
	           !print*, 'Read lonP: ', GD%lonP
	        case ('latp')
	           read(buffer, *, iostat=ios) GD%latP
	           !print*, 'Read latP: ', GD%latP
!********************************************************************
	        case default
	           !print*, 'Skipping invalid label at line', line
	        end select
        end if
     end if
  end do

! Fix for the case when np is absent
if (GD%np == 0) then
	GD%np = GD%nu + GD%nv + GD%nh
end if

close(fh)

end subroutine

!**********************************************************************************************************

subroutine sht_smooth_H(lon, lat, H, P, lon_sht, lat_sht, H_smoothed, topo_dir, topo_file)

!% Given H on a (phi,th) grid, calculates H_smoothed

implicit none

     real(dp), intent(in)		:: lat(:), lon(:), H(:,:)
     character(len=*), intent(in)  :: topo_dir, topo_file
     type(params), intent(in)	:: P

     real(dp), allocatable	:: H_smoothed(:,:), smooth(:)
     integer			    :: nlon, nlat

     real(dp), allocatable	:: lon_sht(:), lat_sht(:)
     integer			    :: nlon_sht, nlat_sht
     real(dp)               :: dlon_sht, dlat_sht

     integer            :: ntrunc, istat, j, deln
     character(len = *), parameter :: gridtype = 'REG'	! (regularly spaced grid)
!     real    ::      T1, T2
	 logical, parameter :: save_coeffs = .true. ! choose between a slower method that saves SHT coeffs for the future usage and quicker method that doesn't

! shortcuts
nlon = size(lon)
nlat = size(lat)

if ( (nlon /= size(H, 1)).or.(nlat /= size(H, 2)) ) then
	print *, 'Dimensions of lon, lat and H do not agree'
	stop
end if

!***************************************************
! On a regular grid spanning -pi/2<th<pi/2 and 0<ph<2*pi
! the condition (ntrunc <= nlat-1) is equiv. to (2*ntrunc <= nlon)
	ntrunc = nlon/2

!***************************************************
! Construct a regular spaced lon-lat grid for the SHT routine
! as the basis take the lon partition of [0, 2*pi) interval (lon_sht = lon)
! Lat coord lat_sht must span [-pi/2, pi/2] and must start at the north pole

nlon_sht = 2*ntrunc
dlon_sht = 2*pi/nlon_sht
nlat_sht = nlon_sht/2 + 1
dlat_sht = pi/(nlat_sht-1)

if (dlon_sht /= dlat_sht) then
	print *, 'dlon_sht /= dlat_sht in SHT subroutine calc_hsal_gridded'
	stop
end if
!print *, "hi1", nlon, nlon_sht, nlat, nlat_sht

! initialises the sht grid structure
    allocate(lon_sht(nlon_sht), stat = istat)
    lon_sht = (/ ( (2*pi/nlon_sht*j), j=0,nlon_sht-1 ) /)	! the grid for the SHT routine spans longitude points
! Th must be the colatitude needed by the SHT routine (increases from the North pole southwards)
    allocate(lat_sht(nlat_sht), stat = istat)
    lat_sht = (/ ( (pi/2 - pi/(nlat_sht-1)*j), j=0,nlat_sht-1 ) /)

!print *, "hi2", lon_sht(1), lon_sht(nlon_sht), lon(1), lon(nlon_sht)
!print *, "hi3", lat_sht(1), lat_sht(nlat_sht), lat(1), lat(nlat_sht)

! Check if the original array satisfies the requirements or if interpolation has to be applied
if ( (lon_sht(1) == lon(1)).and.(lon_sht(nlon_sht) == lon(nlon_sht)) ) then
	if ( (lat_sht(1) == lat(1)).and.(lat_sht(nlat_sht) == lat(nlat_sht)) ) then
		allocate(H_smoothed(nlon_sht, nlat_sht), stat = istat)
		H_smoothed(:,:) = H(1:nlon_sht, 1:nlat_sht)
	elseif ( (lat_sht(1) == lat(nlat_sht)).and.(lat_sht(nlat_sht) == lat(1)) ) then
		allocate(H_smoothed(nlon_sht, nlat_sht), stat = istat)
		H_smoothed(:,:) = H(1:nlon_sht, nlat_sht:1:-1)
		print *, "the right thing"
	else
!		print *, "hi4"
		call bilinear_2d(nlon, nlat, lon, lat, H, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed, P%lib)
	endif
else
! otherwise interpolate H to the grid
    call bilinear_2d(nlon, nlat, lon, lat, H, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed, P%lib)
endif

!	call save_topo(topo_dir // topo_file, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed)

!Smoothing as a function of total wavenumber
allocate(smooth(nlat_sht), stat = istat)
smooth = 0.

if ((P%sht_smooth_H == 0).or.(P%sht_smooth_H >= nlat_sht)) then
		deln = nlat_sht
else
		deln = P%sht_smooth_H + 1
endif

if (P%smooth_type==1) then
	do j = 1, deln
		smooth(j) = 1. ! straightforward truncation smoothing;
	enddo
elseif (P%smooth_type==2) then
	do j = 1, nlat_sht
		smooth(j) = exp(-(float(j-1)/deln)**2)!Gaussian spectral smoothing, where deln is the e-folding width of the smoother in total wavenumber space.
	enddo
elseif (P%smooth_type==3) then
	do j = 1, deln
		smooth(j) = sinc(float(j-1)/deln) !	sinc filter in the spectral space mitigates the Gibbs phenomenon
	enddo
else

	print *, "Unexpected SHT smoothing type. Error in function generate_global_grid."
	stop
endif

! Do the main calculation to find smoothed H on the grid
!	print *, maxval(H_smoothed),minval(H_smoothed)
if (save_coeffs) then
	call specsmooth(H_smoothed,smooth,gridtype, (P%messages > 0), topo_dir//'sht_coeff_'//topo_file) ! spec_coeff_file
else
	call sht_project(H_smoothed,P%sht_smooth_H)
endif
!	print *, maxval(H_smoothed),minval(H_smoothed)

! Clean up
	call cleanup(gridtype)

if (P%smooth_type>=1) then
   	call save_topo(topo_dir//'smoothed_t'//tostring(P%smooth_type)//'_trunc'//tostring(deln-1)// '_'// &
   	               topo_file, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed)
endif

end subroutine sht_smooth_H

!********************************************************************
!-------------------- OLD VERSION, OBSOLETE -------------------------
!********************************************************************
    subroutine interp_H_grid_old(numLons, numLats, xValues, yValues, zValues, &
        numLonsI, numLatsI, xValuesI, yValuesI, zValuesI)


     integer ::      status
     integer ::      i, j
     integer, allocatable :: indi(:), indj(:)
     real(dp)  ::      dx, dy!, dxI, dyI

     integer     ::     numLons, numLats
     real(dp), allocatable, intent(in) :: xValues(:), yValues(:)
     integer, allocatable, intent(in) :: zValues(:, :)

     real(dp)     ::     weight1, weight2
     real(dp)     ::     ztmp1, ztmp2
     integer     ::     numLonsI, numLatsI
     real(dp), allocatable, intent(in) :: xValuesI(:), yValuesI(:)
     integer, allocatable :: zValuesI(:, :)


write(*, '("Interpolating topo to ", f5.2, " min resolution global grid: ")', advance = 'no') 360*60/real(numLonsI)

If (numLons <= numLonsI .or. numLats <= numLatsI .or. numLonsI <= 0 .or. numLatsI <= 0 ) then
     Write (*, *) 'Bad Dimenisons...subroutine interp_H_grid 1'
     print *, numLons, numLonsI, " and ", numLats, numLatsI
     stop
end if

If (xValues(1) >= xValues(numLons) .or. yValues(1) >= yValues(numLats) ) then
     Write (*, *) 'Bad Inputs...subroutine interp_H_grid 2'
     stop
end if

if (allocated(ZValuesI)) deallocate(ZValuesI)
    allocate(ZValuesI(numLonsI, numLatsI), stat = status)
!    if (allocated(XValuesI)) deallocate(XValuesI)
!    allocate(XValuesI(numLonsI), stat = statusx)
!    if (allocated(YValuesI)) deallocate(YValuesI)
!    allocate(YValuesI(numLatsI), stat = statusy)

!    if (statusx /= 0 .or. statusy /= 0 .or. statusz /= 0) then
    if (status /= 0) then
         print *, "[Topo] Can't allocate the memory."
        stop
    end if

    dx = (xValues(numLons) - xValues(1))/(numLons - 1) !the topography is periodical, i.e. map from -180 to 180
    dy = (yValues(numLats) - yValues(1))/(numLats - 1) !and data for xValues(1) and xValues(numLons) is the same

!    dxI = (xValues(numLons) - xValues(1))/numLonsI !grid points are NOT on the -180~180 boundary
!    dyI = (yValues(numLats) - yValues(1))/(numLatsI - 1)

!XValuesI = (/ ( (xValues(1) + dxI*(i-1)), i=1,numLonsI ) /)
!yValuesI = (/ ( (yValues(1) + dyI*(j-1)), j=1,numLatsI ) /)

!Index to know which points of the original data use for interpolation
allocate(indi(numLonsI), stat = status)
allocate(indj(numLatsI), stat = status)

do i = 1, numLonsI
indi(i) = (xValuesI(i)-xValues(1))/dx + 1
enddo
do j = 1, numLatsI
indj(j) = (yValuesI(j)-yValues(1))/dy + 1
enddo

    do i = 1, numLonsI
         do j = 1, numLatsI
         !     First horizontal averaging
          weight1 = (xValuesI(i) - xValues(indi(i)))/dx
          weight2 = (xValues(indi(i)+1) - xValuesI(i))/dx
          if ( (weight1<0 - eps) .or. (weight1>1 + eps) .or. (weight2<0 - eps) .or. (weight2>1 + eps) )  then
                call handle_av_err('Incorrect averaging 1', real(weight1,kind=dp), real(weight2,kind=dp));
                stop
          end if

          ztmp1 = ZValues(indi(i) + 1, indj(j))*weight1 + &
                              ZValues(indi(i), indj(j))*weight2
          ztmp2 = ZValues(indi(i) + 1, indj(j)+1)*weight1 + &
                              ZValues(indi(i), indj(j)+1)*weight2

         !     Then vertical
          weight1 = (yValuesI(j) - yValues(indj(j)))/dy
          weight2 = (yValues(indj(j)+1) - yValuesI(j))/dy
          if (weight1<0 - eps .or. weight1>1 + eps .or. weight2<0 - eps .or. weight2>1 + eps )  then
                call handle_av_err('Incorrect averaging 2', real(weight1,kind=dp), real(weight2,kind=dp));
                stop
          end if

          ZValuesI(i, j) = ztmp2*weight1 + ztmp1*weight2
         enddo
    enddo

        print *, "done."

end subroutine
!***********************************************************************************************************

end module generate_grid



