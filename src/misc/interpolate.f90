module interpolate

     use precisions, only: sp, dp, wp, cwp
     use my_trigs
     use my_sparse
     use save_load
!     use my_sparse_aggregate

!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
     interface mv_average
          module procedure mv_avrg
          module procedure mv_avrg_cmplx
     end interface
!------------------- BILINEAR INTERPOLATION (GRID TO GRID)-------------------------------
     interface bilinear_2d
          module procedure bilinear_grid2grid_4
          module procedure bilinear_grid2grid_8
          module procedure bilinear_grid2grid_cmplx
     end interface

     interface bilinear_2d_mats ! careful. not general, was tweaked for a very specific purpose
		  module procedure bilinear_grid2grid_mats
     end interface

     interface add_border_for_interp
          module procedure add_border_for_interp_
          module procedure add_border_for_interp_cmplx
     end interface

     contains

!====================================================================
!---------- SMOOTHING ROUTINES: mv_average & calc_sm_beta -----------
!====================================================================
subroutine mv_avrg_cmplx(dit, mask, metrics, nf, coor, lib)

!% Smooths a matrix through the weighted moving average method.

implicit none

     complex(cwp),allocatable,intent(inout)	:: dit(:,:)
     integer, intent(in)				    :: nf, coor
     character, intent(in)				    :: lib
     integer	    :: nph, nth

	! averaging
     complex(cwp), allocatable	:: tmp1(:,:)!, tmp2(:,:)
     real(wp), allocatable	:: tmp3(:,:)!, tmp4(:,:)

     real(wp), allocatable	:: metrics(:, :)!, th_hg(:), ta_ug(:), ph_hg(:)
     integer				:: mask(:,:)

     integer            :: j, c, istat, ntrunc, nnz_ph, ncurrent, cth
     integer, allocatable :: nf_ph(:), nf_th(:), ind_th_l(:), ind_th_r(:), nnz_th(:)
     type (triplet_int)	:: frontmat_ph!, backmat_th

!     character(len = *) :: dir_grid, dir_cols

!***************************************************
!% dims
nph = size(dit, 1)
nth = size(dit, 2)

! NOTE: the dims of the averaging box depend on the latitude;
! Mercator (coor = 1): factor in both directions by (1/metrics)^.5
! Lat-lon  (coor = 2): factor in ph-direction by (1/metrics)

! Count the number of elements in the kernel matrices
allocate(nf_ph(nth), nf_th(nth), stat=istat)

	!	for averaging in ph
		if (coor == 1) then
			nf_ph(:)=floor( nf / metrics(1, :)**0.5 )
		elseif (coor == 2) then
			nf_ph(:)=floor( nf / metrics(1, :) )
		endif

	!	for averaging in th
		if (coor == 1) then
			nf_th(:)=floor( nf / metrics(1, :)**0.5 )
		elseif (coor == 2) then
			nf_th(:)=nf
		endif
!		print *, 'nf_ph', minval(nf_ph), maxval(nf_ph)
!		print *, 'nf_th', minval(nf_th), maxval(nf_th)

		! Sparse mat for averaging in th does not change within the th-loop. Initialize just once

		!	first count the number of elements
		allocate(ind_th_l(nth), ind_th_r(nth), nnz_th(nth), stat=istat)
			do c = 1,nth
				ind_th_l(c) = max(1,c-nf_th(c))
				ind_th_r(c) = min(c+nf_th(c),nth)
				nnz_th(c) = ind_th_r(c) - ind_th_l(c) + 1
			end do

!			call init_sparse_0(backmat_th, nth, nth, sum(nnz_th))
!			ncurrent = 1
!			do c = 1,nth
!				backmat_th%indi(ind_th_l(c) : ind_th_r(c)) = (/( j, j=ind_th_l(c), ind_th_r(c) )/)
!				backmat_th%indj(ind_th_l(c) : ind_th_r(c)) = c
!				backmat_th%vals(ind_th_l(c) : ind_th_r(c)) = 1
!				ncurrent = ncurrent + nnz_th(c)
!			end do
!***************************************************************************
! BC the dims of the averaging box depend on the latitude, do a loop in th
!***************************************************************************
		do cth = 1,nth

		!	averaging over 2*nf+1 cells (extra nf in ph and th directions)
			nnz_ph=nph*(2*nf_ph(cth)+1)
    		call init_sparse_0(frontmat_ph, nph, nph, nnz_ph)

			ncurrent = 1
			do c = 1,nph
				frontmat_ph%indi(ncurrent : ncurrent + 2*nf_ph(cth)) = c
				frontmat_ph%indj(ncurrent : ncurrent + 2*nf_ph(cth)) = 1 + modulo( (/(j, j=c-nf_ph(cth),c+nf_ph(cth))/)-1, nph)
				frontmat_ph%vals(ncurrent : ncurrent + 2*nf_ph(cth)) = 1

				ncurrent = ncurrent + 2*nf_ph(cth)+1
			end do

	!***************************************************
			allocate(tmp1(nph, nnz_th(cth)), stat = istat) ! complex valued
!			allocate(tmp2(nph, 1), stat = istat)
			call coo_mat_mul(frontmat_ph, nph, nnz_th(cth), mask(:, ind_th_l(cth) : ind_th_r(cth))*&
							metrics(:, ind_th_l(cth) : ind_th_r(cth))*dit(:, ind_th_l(cth) : ind_th_r(cth)), lib, tmp1)
!			call mat_coo_mul(nph, nth, tmp1,backmat_th, lib, tmp2)
!			dit = tmp2
			dit(:, cth) = sum(tmp1,DIM=2)
			deallocate(tmp1)
!			deallocate(tmp2)

			allocate(tmp3(nph, nnz_th(cth)), stat = istat) ! real valued
!			allocate(tmp4(nph, 1), stat = istat)
			call coo_mat_mul(frontmat_ph, nph, nnz_th(cth), mask(:, ind_th_l(cth) : ind_th_r(cth))*&
							 metrics(:, ind_th_l(cth) : ind_th_r(cth)), lib, tmp3)
!			call mat_coo_mul(nph, nth, tmp3,backmat_th, lib, tmp4)
		!	normalize and apply mask
			where (mask(:,cth) .ne. 0)
!				dit = dit/tmp4
				dit(:, cth) = dit(:, cth)/sum(tmp3,DIM=2)
			elsewhere
				dit(:, cth) = 0
			end where
			deallocate(tmp3)
!			deallocate(tmp4)
!***************************************************************************
		enddo

end subroutine mv_avrg_cmplx

subroutine mv_avrg(dit, mask, metrics, nf, coor, lib)

!% Smooths a matrix through the weighted moving average method.

implicit none

     real(wp),allocatable,intent(inout)	:: dit(:,:)
     integer, intent(in)				    :: nf, coor
     character, intent(in)				    :: lib
     integer	    :: nph, nth

	! averaging
     real(wp), allocatable	:: tmp1(:,:)!, tmp2(:,:)
     real(wp), allocatable	:: tmp3(:,:)!, tmp4(:,:)

     real(wp), allocatable	:: metrics(:, :)!, th_hg(:), ta_ug(:), ph_hg(:)
     integer				:: mask(:,:)

     integer            :: j, c, istat, ntrunc, nnz_ph, ncurrent, cth
     integer, allocatable :: nf_ph(:), nf_th(:), ind_th_l(:), ind_th_r(:), nnz_th(:)
     type (triplet_int)	:: frontmat_ph!, backmat_th

!     character(len = *) :: dir_grid, dir_cols

!***************************************************
!% dims
nph = size(dit, 1)
nth = size(dit, 2)

! NOTE: the dims of the averaging box depend on the latitude;
! Mercator (coor = 1): factor in both directions by (1/metrics)^.5
! Lat-lon  (coor = 2): factor in ph-direction by (1/metrics)

! Count the number of elements in the kernel matrices
allocate(nf_ph(nth), nf_th(nth), stat=istat)

	!	for averaging in ph
		if (coor == 1) then
			nf_ph(:)=floor( nf / metrics(1, :)**0.5 )
		elseif (coor == 2) then
			nf_ph(:)=floor( nf / metrics(1, :) )
		endif

	!	for averaging in th
		if (coor == 1) then
			nf_th(:)=floor( nf / metrics(1, :)**0.5 )
		elseif (coor == 2) then
			nf_th(:)=nf
		endif
!		print *, 'nf_ph', minval(nf_ph), maxval(nf_ph)
!		print *, 'nf_th', minval(nf_th), maxval(nf_th)

		! Sparse mat for averaging in th does not change within the th-loop. Initialize just once

		!	first count the number of elements
		allocate(ind_th_l(nth), ind_th_r(nth), nnz_th(nth), stat=istat)
			do c = 1,nth
				ind_th_l(c) = max(1,c-nf_th(c))
				ind_th_r(c) = min(c+nf_th(c),nth)
				nnz_th(c) = ind_th_r(c) - ind_th_l(c) + 1
			end do

!			call init_sparse_0(backmat_th, nth, nth, sum(nnz_th))
!			ncurrent = 1
!			do c = 1,nth
!				backmat_th%indi(ind_th_l(c) : ind_th_r(c)) = (/( j, j=ind_th_l(c), ind_th_r(c) )/)
!				backmat_th%indj(ind_th_l(c) : ind_th_r(c)) = c
!				backmat_th%vals(ind_th_l(c) : ind_th_r(c)) = 1
!				ncurrent = ncurrent + nnz_th(c)
!			end do
!***************************************************************************
! BC the dims of the averaging box depend on the latitude, do a loop in th
!***************************************************************************
		do cth = 1,nth

		!	averaging over 2*nf+1 cells (extra nf in ph and th directions)
			nnz_ph=nph*(2*nf_ph(cth)+1)
    		call init_sparse_0(frontmat_ph, nph, nph, nnz_ph)

			ncurrent = 1
			do c = 1,nph
				frontmat_ph%indi(ncurrent : ncurrent + 2*nf_ph(cth)) = c
				frontmat_ph%indj(ncurrent : ncurrent + 2*nf_ph(cth)) = 1 + modulo( (/(j, j=c-nf_ph(cth),c+nf_ph(cth))/)-1, nph)
				frontmat_ph%vals(ncurrent : ncurrent + 2*nf_ph(cth)) = 1

				ncurrent = ncurrent + 2*nf_ph(cth)+1
			end do

	!***************************************************
			allocate(tmp1(nph, nnz_th(cth)), stat = istat) ! real valued
!			allocate(tmp2(nph, 1), stat = istat)
			call coo_mat_mul(frontmat_ph, nph, nnz_th(cth), mask(:, ind_th_l(cth) : ind_th_r(cth))*&
							metrics(:, ind_th_l(cth) : ind_th_r(cth))*dit(:, ind_th_l(cth) : ind_th_r(cth)), lib, tmp1)
!			call mat_coo_mul(nph, nth, tmp1,backmat_th, lib, tmp2)
!			dit = tmp2
			dit(:, cth) = sum(tmp1,DIM=2)
			deallocate(tmp1)
!			deallocate(tmp2)

			allocate(tmp3(nph, nnz_th(cth)), stat = istat) ! real valued
!			allocate(tmp4(nph, 1), stat = istat)
			call coo_mat_mul(frontmat_ph, nph, nnz_th(cth), mask(:, ind_th_l(cth) : ind_th_r(cth))*&
							 metrics(:, ind_th_l(cth) : ind_th_r(cth)), lib, tmp3)
!			call mat_coo_mul(nph, nth, tmp3,backmat_th, lib, tmp4)
		!	normalize and apply mask
			where (mask(:,cth) .ne. 0)
!				dit = dit/tmp4
				dit(:, cth) = dit(:, cth)/sum(tmp3,DIM=2)
			elsewhere
				dit(:, cth) = 0
			end where
			deallocate(tmp3)
!			deallocate(tmp4)
!***************************************************************************
		enddo

end subroutine mv_avrg
!**********************************************************************************************************

subroutine calc_sm_beta(dit, h, mask, metrics, nf, beta, lib)

!% Smooths a matrix through the weighted moving average method.

implicit none

     complex(cwp),allocatable,intent(in) 	:: dit(:,:), h(:,:)
     integer, intent(in)				    :: nf
     character, intent(in)				    :: lib
     integer							    :: nf_local
     complex(cwp), allocatable				:: beta(:,:)
     complex(cwp)							:: beta0
     integer	    :: nph, nth

	! averaging
     complex(cwp), allocatable	:: tmp1(:,:), tmp2(:,:)!, beta_tmp(:,:)
     real(wp), allocatable	:: tmp3(:,:), tmp4(:,:)

     real(wp), allocatable	:: metrics(:, :)!, th_hg(:), ta_ug(:), ph_hg(:)
     integer				:: mask(:,:)

     integer            :: j, c, istat, ntrunc, nnz, ncurrent
     type (triplet_int)	:: frontmat_ph, backmat_th
     complex(cwp), allocatable	:: beta_u(:,:), beta_v(:,:)

!     character(len = *) :: dir_grid, dir_cols

!***************************************************
!% dims
nph = size(h, 1)
nth = size(h, 2)

	!***************************************************
	allocate(beta(nph, nth), stat = istat)
		beta0 = sum( mask*metrics*dit*conjg(h) ) / sum( mask*metrics*h*conjg(h) )
		beta = beta0

		!	averaging over 2*nf+1 cells (extra nf in ph and th directions)
			nnz=nph*(2*nf+1)
    		call init_sparse_0(frontmat_ph, nph, nph, nnz)

			ncurrent = 1
			do c = 1,nph
				frontmat_ph%indi(ncurrent : ncurrent + 2*nf) = c
				frontmat_ph%indj(ncurrent : ncurrent + 2*nf) = 1 + modulo( (/(j, j=c-nf,c+nf)/)-1, nph)
				frontmat_ph%vals(ncurrent : ncurrent + 2*nf) = 1

				ncurrent = ncurrent + 2*nf+1
			end do

		!	averaging in th
		!	first count the number of elements
			nnz = 0
			do c = 1,nth
				! Mercator not implemented
!				nf_local=ceiling( nf /cos(th_hg(c)) ) ! 1/cos(th_ug(cta)) bc of Mercator coords
				nf_local=nf
				nnz = nnz + min(c+nf_local,nth) - max(1,c-nf_local) + 1
			end do
			call init_sparse_0(backmat_th, nth, nth, nnz)

			ncurrent = 1
			do c = 1,nth
				! Mercator not implemented
!				nf_local=ceiling( nf /cos(th_hg(c)) ) ! 1/cos(th_ug(cta)) bc of Mercator coords
				nf_local=nf
				nnz = min(c+nf_local,nth) - max(1,c-nf_local) + 1

				backmat_th%indi(ncurrent : ncurrent + (nnz-1)) = (/(j, j=max(1,c-nf_local),min(c+nf_local,nth))/)
				backmat_th%indj(ncurrent : ncurrent + (nnz-1)) = c
				backmat_th%vals(ncurrent : ncurrent + (nnz-1)) = 1

				ncurrent = ncurrent + nnz
			end do

!			allocate(beta(nph, nth), stat = istat)

	!***************************************************
			allocate(tmp1(nph, nth), tmp2(nph, nth), stat = istat) ! complex valued
			call coo_mat_mul(frontmat_ph, nph, nth, metrics*dit*conjg(h), lib, tmp1)
			call mat_coo_mul(nph, nth, tmp1,backmat_th, lib, tmp2)
			beta = tmp2
			deallocate(tmp1, tmp2)

			allocate(tmp3(nph, nth),tmp4(nph, nth), stat = istat) ! real valued
			call coo_mat_mul(frontmat_ph, nph, nth, metrics*(abs(h)**2), lib, tmp3)
			call mat_coo_mul(nph, nth, tmp3,backmat_th, lib, tmp4)
			where (mask .ne. 0)
				beta = beta / tmp4
			elsewhere
				beta = 0
			end where

!				print *, "divide"
!				print *, count(isnan(real(beta))), count((mask .ne. 0).and.(tmp4==0))

			deallocate(tmp3, tmp4)


end subroutine calc_sm_beta


!====================================================================
!---------- INTERPOLATEION ON A SPHERE (RANDOM TO GRID) -----------
!====================================================================
!subroutine save_sol4plot_cmplx(nph, nth, nvar, varp, var, sol4plot_name)
!
!    implicit none
!
!     integer, intent(in) 		:: nph, nth, nvar
!     integer, intent(in) 		:: varp(:, :)
!     complex(cwp), intent(in)	:: var(:)
!
!    character(len=*) :: sol4plot_name
!
!      OPEN(10,status='unknown', file=sol4plot_name, form='unformatted', action='write', access='stream')
!
!      write(10) nph, nth, nvar
!      write(10) varp(:,1)
!      write(10) varp(:,2)
!      write(10) abs(var)
!      CLOSE(10)
!
!
!end subroutine save_sol4plot_cmplx
!====================================================================
!------- BILINEAR INTERPOLATION COMPLEX (GRID TO GRID) --------------
!====================================================================
subroutine bilinear_grid2grid_cmplx(nph, nth, ph, th, f, nph1, nth1, ph1, th1, f1, lib)!, dir_grid)
! (use simple bilinear interpolation)

!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph < 2*pi,
!% interpolates to the grid {phi1, theta1}. The
!% output grid is specified by the two vectors
!% (th1,ph1).
!% Assumes that th, ph, th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none
!character(len=*) :: dir_grid

integer, intent(in)  :: nph, nth, nph1, nth1
complex(cwp),intent(in)  :: f(nph, nth) ! old data
real(wp),intent(in)  :: ph(nph), th(nth), ph1(nph1), th1(nth1)
character, intent(in)	:: lib

complex(wp),allocatable :: f1(:,:) ! interpolated data

type (triplet)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
complex(cwp), allocatable:: vec_tmp(:), mat_tmp(:, :)

integer:: c, loc, istatus

! Prepare interpolation matrices
!print *, nph, nth, nph1, nth1

!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.
call init_sparse_0(im_th, nth1, nth, 2*nth1)
do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!print *, ""

!% given a grid ph with 0 <= ph < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi
call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate(nph, ph, ph1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
!  print *, im_ph%vals(2*c-1),im_ph%vals(2*c)
end do
! first average in ph
    allocate(vec_tmp(nph1), stat = istatus)
    allocate(mat_tmp(nth, nph1), stat = istatus)

do c = 1, nth
	call coo_vec_mul(im_ph,f(:, c),lib, vec_tmp)
	mat_tmp(c, :) = vec_tmp
end do

! then average in th
    if (allocated(vec_tmp)) deallocate(vec_tmp)
    allocate(vec_tmp(nth1), stat = istatus)

    if (allocated(f1)) deallocate(f1)
    allocate(f1(nph1, nth1), stat = istatus)

do c = 1, nph1
	call coo_vec_mul(im_th,mat_tmp(:, c),lib, vec_tmp)
	f1(c, :) = vec_tmp
end do

!call save_cmat(im_ph, dir_grid // 'im_ph.dat')
!call save_cmat(im_th, dir_grid // 'im_th.dat')

end subroutine bilinear_grid2grid_cmplx

!====================================================================
!---------- BILINEAR INTERPOLATION (GRID TO GRID) ------------------
!====================================================================
subroutine bilinear_grid2grid_4(nph, nth, ph, th, f, nph1, nth1, ph1, th1, f1,lib)!, dir_grid)
! (use simple bilinear interpolation)

!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph < 2*pi,
!% interpolates to the grid {phi1, theta1}. The
!% output grid is specified by the two vectors
!% (th1,ph1).
!% Assumes that th, ph, th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none
!character(len=*) :: dir_grid

integer, intent(in)  :: nph, nth, nph1, nth1
real(sp),intent(in)  :: f(nph, nth) ! old data
real(sp),intent(in)  :: ph(nph), th(nth), ph1(nph1), th1(nth1)
character, intent(in)	:: lib

real(sp),allocatable :: f1(:,:) ! interpolated data

type (triplet_sp)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
real(sp), allocatable:: vec_tmp(:), mat_tmp(:, :)

integer:: c, loc, istatus

! Prepare interpolation matrices
!print *, nph, nth, nph1, nth1
!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.

call init_sparse_0(im_th, nth1, nth, 2*nth1)

do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate_4(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!% given a grid ph with 0 <= ph < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi
call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate_4(nph, ph, ph1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
!  print *, im_ph%vals(2*c-1),im_ph%vals(2*c)
end do
! first average in ph
    allocate(vec_tmp(nph1), stat = istatus)
    allocate(mat_tmp(nth, nph1), stat = istatus)

do c = 1, nth
	call coo_vec_mul(im_ph,f(:, c),lib, vec_tmp)
	mat_tmp(c, :) = vec_tmp
end do

! then average in th
    if (allocated(vec_tmp)) deallocate(vec_tmp)
    allocate(vec_tmp(nth1), stat = istatus)

    if (allocated(f1)) deallocate(f1)
    allocate(f1(nph1, nth1), stat = istatus)

do c = 1, nph1
	call coo_vec_mul(im_th,mat_tmp(:, c),lib, vec_tmp)
	f1(c, :) = vec_tmp
end do

!call save_cmat(im_ph, dir_grid // 'im_ph.dat')
!call save_cmat(im_th, dir_grid // 'im_th.dat')

end subroutine bilinear_grid2grid_4

subroutine bilinear_grid2grid_8(nph, nth, ph, th, f, nph1, nth1, ph1, th1, f1, lib)!, dir_grid)
! (use simple bilinear interpolation)

!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph < 2*pi,
!% interpolates to the grid {phi1, theta1}. The
!% output grid is specified by the two vectors
!% (th1,ph1).
!% Assumes that th, ph, th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none
!character(len=*) :: dir_grid

integer, intent(in)  :: nph, nth, nph1, nth1
real(dp),intent(in)  :: f(nph, nth) ! old data
real(dp),intent(in)  :: ph(nph), th(nth), ph1(nph1), th1(nth1)
character, intent(in)	:: lib

real(dp),allocatable :: f1(:,:) ! interpolated data

type (triplet_dp)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
real(dp), allocatable:: vec_tmp(:), mat_tmp(:, :)

integer:: c, loc, istatus

! Prepare interpolation matrices
!print *, nph, nth, nph1, nth1
!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.
call init_sparse_0(im_th, nth1, nth, 2*nth1)

do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate_8(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!% given a grid ph with 0 <= ph < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi
call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate_8(nph, ph, ph1(c), loc)
!	  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
!  print *, im_ph%vals(2*c-1),im_ph%vals(2*c)
end do
! first average in ph
    allocate(vec_tmp(nph1), stat = istatus)
    allocate(mat_tmp(nth, nph1), stat = istatus)

do c = 1, nth
	call coo_vec_mul(im_ph,f(:, c),lib, vec_tmp)
	mat_tmp(c, :) = vec_tmp
end do

! then average in th
    if (allocated(vec_tmp)) deallocate(vec_tmp)
    allocate(vec_tmp(nth1), stat = istatus)

    if (allocated(f1)) deallocate(f1)
    allocate(f1(nph1, nth1), stat = istatus)

do c = 1, nph1
	call coo_vec_mul(im_th,mat_tmp(:, c),lib, vec_tmp)
	f1(c, :) = vec_tmp
end do

!call save_cmat(im_ph, dir_grid // 'im_ph.dat')
!call save_cmat(im_th, dir_grid // 'im_th.dat')

end subroutine bilinear_grid2grid_8
!**********************************************************************************************************
subroutine bilinear_grid2grid_mats(nph, nth, ph, th, nph1, nth1, ph1, th1, im_ph, im_th)
! (return mats for bilinear interpolation)
!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph (mod 2*pi) < 2*pi,

!% Assumes that th, ph (mod 2*pi), th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none

integer, intent(in)  :: nph, nth, nph1, nth1
real(wp),intent(in)  :: th(nth), ph1(nph1), th1(nth1)
real(wp)			 :: ph(nph)

type (triplet)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
real(wp), allocatable:: dph_tmp(:)

integer:: c, loc, istatus

! Prepare interpolation matrices
!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.
call init_sparse_0(im_th, nth1, nth, 2*nth1)

do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!% given a grid ph with 0 <= ph (mod 2*pi) < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi

! Take care of situations where ph is monotonic but only after adjusting to (mod 2*pi)
! I.e., neighbours used for interpolations were sticking out (either the first or the last one)
! WARNING: this is very much ad-hoc
! Check monotonicity:
allocate(dph_tmp(nph), stat=istatus)
dph_tmp = cshift(ph, shift=1)-ph
if (any(dph_tmp(1:nph-1)<0)) then
	if (dph_tmp(1)<0) then
		ph(1) = ph(1)-2*pi
	elseif (dph_tmp(nph-1)<0) then
		ph(nph) = ph(nph)+2*pi
	else
		print *, "Error in bilinear_grid2grid_mats:"
		print *, " ph  ", ph
		print *, " dph ", dph_tmp
		stop
	endif
	! check again that this helped
	dph_tmp = cshift(ph, shift=1)-ph
	if (any(dph_tmp(1:nph-1)<0)) then
		print *, "Checked again. Error in bilinear_grid2grid_mats:"
		print *, " ph  ", ph
		print *, " dph ", dph_tmp
		stop
	endif
endif

call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate(nph, ph, ph1(c), loc)
!	  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
end do

end subroutine bilinear_grid2grid_mats
!**********************************************************************************************************
subroutine locate(n, xx, x, j)
!% given an array xx(1:n), and a value x, returns a value
!% j such that x is between xx(j) and xx(j+1). xx(1:n) must
!% be monotonic, either increasing or decreasing. j=0 or j=n
!% is returned to indicate that x is out of range.
    implicit none

integer, intent(in)  :: n
real(wp),intent(in)  :: xx(n), x

integer :: j, jl, ju, jm

! First consider the out of range cases
if ( (xx(n)>=xx(1)).and.(x <= xx(1)) .or. (xx(n)<=xx(1)).and.(x >= xx(1)) ) then
  j=0
elseif ( (xx(n)>=xx(1)).and.(x >= xx(n)) .or. (xx(n)<=xx(1)).and.(x <= xx(n)) ) then
  j=n
else

	jl=0
	ju=n+1

	do while (ju-jl > 1)
	  jm =(ju+jl)/2 ! implicit floor rounding
	  if ( (xx(n)>=xx(1)) .eqv. (x>=xx(jm)) ) then
	     jl=jm
	  else
	    ju=jm
	  end if
	end do

	j = jl

end if
!
!if (j == 0) then
!print *, x
!endif

end subroutine locate

subroutine locate_4(n, xx, x, j)

    implicit none

integer, intent(in)  :: n
real(sp),intent(in)  :: xx(n), x

integer :: j, jl, ju, jm

! First consider the out of range cases
if ( (xx(n)>=xx(1)).and.(x <= xx(1)) .or. (xx(n)<=xx(1)).and.(x >= xx(1)) ) then
  j=0
elseif ( (xx(n)>=xx(1)).and.(x >= xx(n)) .or. (xx(n)<=xx(1)).and.(x <= xx(n)) ) then
  j=n
else

	jl=0
	ju=n+1

	do while (ju-jl > 1)
	  jm =(ju+jl)/2 ! implicit floor rounding
	  if ( (xx(n)>=xx(1)) .eqv. (x>=xx(jm)) ) then
	     jl=jm
	  else
	    ju=jm
	  end if
	end do

	j = jl

end if
!
!if (j == 0) then
!print *, x
!endif

end subroutine locate_4

subroutine locate_8(n, xx, x, j)

    implicit none

integer, intent(in)  :: n
real(dp),intent(in)  :: xx(n), x

integer :: j, jl, ju, jm

! First consider the out of range cases
if ( (xx(n)>=xx(1)).and.(x <= xx(1)) .or. (xx(n)<=xx(1)).and.(x >= xx(1)) ) then
  j=0
elseif ( (xx(n)>=xx(1)).and.(x >= xx(n)) .or. (xx(n)<=xx(1)).and.(x <= xx(n)) ) then
  j=n
else

	jl=0
	ju=n+1

	do while (ju-jl > 1)
	  jm =(ju+jl)/2 ! implicit floor rounding
	  if ( (xx(n)>=xx(1)) .eqv. (x>=xx(jm)) ) then
	     jl=jm
	  else
	    ju=jm
	  end if
	end do

	j = jl

end if
!
!if (j == 0) then
!print *, x
!endif

end subroutine locate_8
!**********************************************************************************************************
subroutine add_border_for_interp_(nph, nth, h, mask)
!% given an array h(:,:), and a mask array of the same dims (that indicates the exterior)
!% extends the values of h onto its boundary gridpoints (coastline)
    implicit none

	integer,intent(in)	 :: nph, nth
	real(wp)  :: h(nph, nth)
	logical,intent(in)   :: mask(nph, nth)
	real(wp)			 :: hborder(nph, nth)

	hborder = 0;
    where (mask .and. cshift(.not.(mask), shift=1, dim=1))	hborder = cshift(h, shift=1, dim=1);
    where (mask .and. cshift(.not.(mask), shift=-1, dim=1).and. (hborder==0))	hborder = cshift(h, shift=-1, dim=1);
    where (mask .and. cshift(.not.(mask), shift=1, dim=2).and. (hborder==0))	hborder = cshift(h, shift=1, dim=2);
	where (mask .and. cshift(.not.(mask), shift=-1, dim=2).and. (hborder==0)) hborder = cshift(h, shift=-1, dim=2);
    h = h + hborder;

end subroutine add_border_for_interp_
!**********************************************************************************************************
subroutine add_border_for_interp_cmplx(nph, nth, h, mask)
!% given an array h(:,:), and a mask array of the same dims (that indicates the exterior)
!% extends the values of h onto its boundary gridpoints (coastline)
    implicit none

	integer,intent(in)	 :: nph, nth
	complex(cwp)		 :: h(nph, nth)
	logical,intent(in)   :: mask(nph, nth)
	complex(cwp)		 :: hborder(nph, nth)

	hborder = 0;
    where (mask .and. cshift(.not.(mask), shift=1, dim=1))	hborder = cshift(h, shift=1, dim=1);
    where (mask .and. cshift(.not.(mask), shift=-1, dim=1).and. (hborder==0))	hborder = cshift(h, shift=-1, dim=1);
    where (mask .and. cshift(.not.(mask), shift=1, dim=2).and. (hborder==0))	hborder = cshift(h, shift=1, dim=2);
	where (mask .and. cshift(.not.(mask), shift=-1, dim=2).and. (hborder==0)) hborder = cshift(h, shift=-1, dim=2);
    h = h + hborder;

end subroutine add_border_for_interp_cmplx
!**********************************************************************************************************
!**********************************************************************************************************

end module interpolate
