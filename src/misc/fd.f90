module fd

     use precisions, only: wp
     use my_trigs
     use my_lapack
     use dispmodule
!     use save_load
! my version of DIVDIF Fortran 90 library
!**********************************************************************************************************
!DIVDIF is a FORTRAN90 library which creates, prints and manipulates divided difference polynomials based on data tabulated at evenly spaced or unevenly spaced argument values.
!Divided difference polynomials are a systematic method of computing polynomial approximations to scattered data. The representations are compact, and may easily be updated with new data, rebased at zero, or analyzed to produce the standard form polynomial, integral or derivative polynomials.
!Other routines are available to convert the divided difference representation to standard polynomial format. This is a natural way to determine the coefficients of the polynomial that interpolates a given set of data, for instance.
!One surprisingly simple but useful routine is available to take a set of roots and compute the divided difference or standard form polynomial that passes through those roots.
!Finally, the Newton-Cotes quadrature formulas can be derived using divided difference methods, so a few routines are given which can compute the weights and abscissas of open or closed rules for an arbitrary number of nodes. !

contains

!*****************************************************************************
subroutine savgol_smooth_2d(nph,nth, raw_data, nx, ny, mx, my, sm_data)

integer, intent(in)	:: nph,nth
real(wp), intent(in):: raw_data(nph,nth)

integer, intent(in)	:: nx, ny ! half of the # of pts in x and y
integer, intent(in)	:: mx, my ! polinomial order in x and y
integer :: nlx,nrx,npx, nly,nry,npy

real(wp), allocatable :: sm_data(:,:)
real(wp), allocatable :: c(:,:)
integer :: cth, cph, j, k, cphk, istat

  nlx=nx; nrx=nx; nly=ny; nry=ny;
  npx = nlx+nrx+1	! # of pts used for interpolation
  npy = nly+nry+1	! the idea is to have twice as many pts as needed to uniquely determine a polinomial of degree

  allocate (sm_data(nph,nth), stat=istat)
  sm_data=raw_data ! initialize, for unsmoothed edges around the (top/bottom) boundary

  allocate (c(-nlx:nrx,-nly:nry), stat=istat)

! calculate Savitzky-Golay filter coefficients
  call savgol_2d(c,npx,nlx,nrx,npy,nly,nry,0,0,mx,my, 1) ! 0,0 -- smooth function, 1 - coupling between 1d coeffs
!  print *,' Savitzky-Golay Filter Coefficients:'
!  call disp(c)

! Apply filter to input data
  do cth=1+nly, nth-nry ! skip (top/bottom) boundary points
  	do cph=1, nlx 		! assume periodicity in ph
  		sm_data(cph, cth) = 0.
  		do k=-nlx,nrx
  			cphk=1+modulo((cph+k)-1,nph)
		    sm_data(cph, cth)=sm_data(cph, cth) + sum( c(k,-nly:nry)*raw_data(cphk, cth-nly:cth+nry) )
	    end do
  	end do
  	do cph=1+nph-nrx, nph	! assume periodicity in ph
  		sm_data(cph, cth) = 0.
  		do k=-nlx,nrx
  			cphk=1+modulo((cph+k)-1,nph)
		    sm_data(cph, cth)=sm_data(cph, cth) + sum( c(k,-nly:nry)*raw_data(cphk, cth-nly:cth+nry) )
	    end do
  	end do
  	do cph=1+nlx, nph-nrx
	    sm_data(cph, cth) = sum( c(-nlx:nrx,-nly:nry)*raw_data(cph-nlx:cph+nrx, cth-nly:cth+nry) )
  	end do
  end do

end subroutine

!*****************************************************************************
subroutine savgol_smooth_1d(n, raw_data, nn, m, sm_data)

integer, intent(in)	:: n
real(wp), intent(in):: raw_data(n)

integer, intent(in)	:: nn ! half of the # of pts
integer, intent(in)	:: m ! polinomial order
integer :: nl,nr,np

real(wp), allocatable :: sm_data(:)
real(wp), allocatable :: c(:)
integer :: j, k, istat

  nl=nn; nr=nn;
  np = nl+nr+1	! # of pts used for interpolation
				! the idea is to have twice as many pts as needed to uniquely determine a polinomial of degree m

  allocate (sm_data(n), stat=istat)
  sm_data=raw_data ! initialize, for unsmoothed edges around the (top/bottom) boundary

  allocate (c(-nl:nr), stat=istat)
! calculate Savitzky-Golay filter coefficients
  call savgol_1d(c,np,nl,nr,0,m) ! 0 -- smooth function, 1 - 1st derivative, etc

!  print *,' Savitzky-Golay Filter Coefficients:'
!  do j=1, np
!  	write(*,*) ind(j), c(j)
!  enddo

! Apply filter to input data
  do j=1+nl, n-nr !skip left/right points that do not exist
    sm_data(j)=0.
    do k=-nl,nr
	    sm_data(j)=sm_data(j)+c(k)*raw_data(j+k)
    end do
  end do

end subroutine

!*****************************************************************************
subroutine savgol_2d(h,npx,nlx,nrx,npy,nly,nry,dx,dy,mx,my, flag_coupling)

implicit none

!--------------------------------------------------------------------------------------------
! Same as savgol_1d but in for 2d data
!--------------------------------------------------------------------------------------------
	real(wp), intent(out) 	:: h(-nlx:nrx,-nly:nry)
	integer, intent(in)		:: dx,mx,nlx,nrx,npx
	integer, intent(in)		:: dy,my,nly,nry,npy
	integer, intent(in)		:: flag_coupling

	real(wp)			 	:: h1(-nlx:nrx), h2(-nly:nry)
	integer :: k1,k2

call savgol_1d(h1,npx,nlx,nrx,dx,mx)
call savgol_1d(h2,npy,nly,nry,dy,my)

if (flag_coupling==0) then

    if ((dx==0).and.(dy==0)) then
        do k1=-nlx,nrx
            do k2=-nly,nry
                h(k1,k2)=h1(k1)/npy+h2(k2)/npx-1/npx/npy;
            enddo
        enddo

    elseif (dx==0) then
        do k1=-nlx,nrx
            do k2=-nly,nry
                h(k1,k2)=h2(k2)/npx;
            enddo
        enddo

    elseif (dy==0) then
        do k1=-nlx,nrx
            do k2=-nly,nry
                h(k1,k2)=h1(k1)/npy;
            enddo
        enddo
    else
        h=0
    endif

else
    do k1=-nlx,nrx
        do k2=-nly,nry
            h(k1,k2)=h1(k1)*h2(k2);
        enddo
    enddo
endif

end subroutine savgol_2d
!*****************************************************************************
subroutine savgol_1d(c,np,nl,nr,ld,m)

implicit none

!--------------------------------------------------------------------------------------------
!USES lapack_solve_general given below.
!Returns in c(1:np), in wrap-around order (see reference) consistent with the argument respns
!in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward
!(past) data points used, while nr is the number of rightward (future) data points, making
!the total number of data points used nl +nr+1. ld is the order of the derivative desired
!(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also
!equal to the highest conserved moment; usual values are m = 2 or m = 4.
!--------------------------------------------------------------------------------------------
	real(wp), intent(out) 	:: c(-nl:nr)
	integer, intent(in)		:: ld,m,nl,np,nr

	integer					:: imj,ipj,j,k,mm, istat

	real(wp), allocatable 	:: a(:,:),b(:)
	real(wp)				:: fac,tmp


  if   ((np.lt.nl+nr+1).or.(nl.lt.0).or.(nr.lt.0).or.&
  		(ld.gt.m).or.(nl+nr.lt.m)) then
		print *, ' Bad args in savgol.'; stop;
  endif

	allocate (a(m+1,m+1),b(m+1), stat=istat)

  do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit.
    tmp=0.
    if(ipj.eq.0) tmp=1.
    do k=1,nr
      tmp=tmp+dfloat(k)**ipj
    end do
    do k=1,nl
      tmp=tmp+dfloat(-k)**ipj
    end do
    mm=min(ipj,2*m-ipj)
    do imj=-mm,mm,2
      a(1+(ipj+imj)/2,1+(ipj-imj)/2)=tmp
    end do
  end do

  do j=1,m+1
    b(j)=0.
  end do
  b(ld+1)=1.      !Right-hand side vector is unit vector, depending on which derivative we want.

   call lapack_solve_general(a,b) ! solve the system of linear eqs

  !Zero the output array (it may be bigger than the number of coefficients).
  c=0.

  do k=-nl,nr                        !Each Savitzky-Golay coefficient is the dot product
    tmp=b(1)                         !of powers of an integer with the inverse matrix row.
    fac=1.
    do mm=1,m
      fac=fac*k
      tmp=tmp+b(mm+1)*fac
    end do
    c(k)=tmp
  end do

end subroutine savgol_1d
!*****************************************************************************

subroutine laplace_relax(data, mask)
!*****************************************************************************80
!
! laplace_relax uses the (external!) boundary values of a connected grid-region
!	 and copies them into the interior. Implemented as a propagating boundary front

  implicit none

  integer, intent(in)		:: mask(:,:)
  real(wp), intent(inout)	:: data(:,:)
  integer, parameter  :: iter_max=360 ! max number of iterations (on each of the multiple grids)
  real(wp), parameter :: tol=1.0d-2
!  character(len=*) :: dir_grid
!	successive over-relaxation params:
  real(wp)			  :: omega, rhoj, rhojsq

  integer, parameter  :: min_m_pts = 360 ! min number of x grid points at the coarsest grid
  integer 			  :: cgrid, step, mgrid, ngrid, n_grids

  integer, allocatable	:: maskgrid(:,:)
  real(wp), allocatable :: datagrid(:,:), datagridnew(:,:)
  logical, allocatable	:: ind_b(:,:)
  integer 				:: lb, rb, bb, ub

  real(wp)		:: errij,error, tmp_val
  integer		:: j, l, dims(2), istat, iter, m, n

  if (maxval( abs(shape(mask)-shape(data)) )/=0) then
  	print *, 'Error in laplace_relax:'
  	print *, 'Data and mask arrays must be of the same size'
  	stop
  endif


  dims = shape(mask)
  m = dims(1); n = dims(2);

  n_grids = log(real(m)/min_m_pts)/log(2.) ! number of grids in our grid tree

! Loop from the coarsest to the finest (original) grid
do cgrid=n_grids,0, -1
  step = 2**cgrid
  mgrid = (m-1)/step + 1
  ngrid = (n-1)/step + 1

  allocate ( datagrid(mgrid,ngrid),maskgrid(mgrid,ngrid), stat=istat )
    do j=1,mgrid
      do l=1,ngrid
		maskgrid(j,l) = mask(1+(j-1)*step,1+(l-1)*step)
      end do
    end do

!print *, "mgrid, ngrid, step, nmask", mgrid, ngrid, step, count(maskgrid==1)

    do j=1,mgrid
      do l=1,ngrid
      	if ( (cgrid==n_grids).or.(maskgrid(j,l)==0)) then
			datagrid(j,l) = data(1+(j-1)*step,1+(l-1)*step)
!		elseif ( (modulo(j-1,2)==1).or.(modulo(l-1,2)==1) ) then
!			datagrid(j,l) = min (datagridnew(1+(j-1)/2,1+(l-1)/2),  &
!								 datagridnew(1+modulo(((j-1)/2+2)-1,(mgrid-1)/2+1),1+(l-1)/2), &
!								 datagridnew(1+(j-1)/2,1+(l-1)/2),  &
!								 datagridnew(1+(j-1)/2, min(2+(l-1)/2,(ngrid-1)/2+1)) ) ! min val of 4 surrounding points, uses periodicity
		elseif ( (modulo(j-1,2)==1).and.(modulo(l-1,2)==0) ) then
			datagrid(j,l) = .5* (datagridnew(1+(j-1)/2,1+(l-1)/2)  +  &
								 datagridnew(1+modulo( ((j-1)/2+2)-1,(mgrid-1)/2+1),1+(l-1)/2)) ! average of 4 surrounding points, uses periodicity
		elseif ( (modulo(j-1,2)==0).and.(modulo(l-1,2)==1) ) then
			datagrid(j,l) = .5* (datagridnew(1+(j-1)/2,1+(l-1)/2)  +  &
								 datagridnew(1+(j-1)/2, min(2+(l-1)/2,(ngrid-1)/2+1)) ) ! average of 4 surrounding points, uses periodicity
		elseif ( (modulo(j-1,2)==1).and.(modulo(l-1,2)==1) ) then
			datagrid(j,l) =.25* (datagridnew(1+(j-1)/2,1+(l-1)/2)  +  &
								 datagridnew(1+modulo(((j-1)/2+2)-1,(mgrid-1)/2+1),1+(l-1)/2) + &
								 datagridnew(1+(j-1)/2,1+(l-1)/2)  +  &
								 datagridnew(1+(j-1)/2, min(2+(l-1)/2,(ngrid-1)/2+1)) ) ! average of 4 surrounding points, uses periodicity
		else
			datagrid(j,l) = datagridnew(1+(j-1)/2,1+(l-1)/2)
		endif
      end do
    end do

! 'land' boundary of the sponge (where issponge==0 and cn==0)
  if (allocated(ind_b)) deallocate(ind_b)
  allocate ( ind_b(mgrid,ngrid), stat=istat )
	ind_b = ( (datagrid == 0) .and. (maskgrid == 0) .and. &
	                           ( (cshift(maskgrid, shift=1, dim=1) == 1) .or. &
	                             (cshift(maskgrid, shift=-1, dim=1) == 1) .or. &
	                             (cshift(maskgrid, shift=1, dim=2) == 1) .or. &
	                             (cshift(maskgrid, shift=-1, dim=2) == 1) ) )

! Modify vals of cn at land boundary to be the min values at the same latitude
    do j=1,mgrid
      do l=1,ngrid
		if (ind_b(j,l) .and. any((maskgrid(:, l) == 0).and.(datagrid(:, l) > 0)) ) then ! additional condition for regions around the poles where they may be no non-sponge ocean
			datagrid(j,l) = minval( datagrid(:, l), mask=( (maskgrid(:, l) == 0).and.(datagrid(:, l) > 0) ) )
		endif
      end do
    end do
!	extra loop for regions around the poles where they may be no non-sponge ocean
	tmp_val = minval(datagrid, mask = (ind_b.and.(datagrid>0)))
!	print *, 'tmpval', tmp_val
	where (ind_b .and. (datagrid<=0)) datagrid =tmp_val

  if (allocated(datagridnew)) deallocate(datagridnew)
  allocate ( datagridnew(mgrid,ngrid), stat=istat )
    do j=1,mgrid
      do l=1,ngrid
		datagridnew(j,l) = datagrid(j,l)
      end do
    end do

!    call save_matrix(datagrid, dir_grid // 'itm/' // 'c2n_hg_cgrid_before.dat')

  iter=0
  error=1

!  rhoj = 1 - ((pi/(2*(m+1))) * (pi/(2*(m+1))))
!  rhojsq = rhoj*rhoj
  omega = 1
!  omega = 1/(1-rhojsq/2)

  do while ( error .gt. tol/step .and. iter .lt. iter_max )
    error=0.0d0

    do j=1,mgrid
      do l=1,ngrid
		if (maskgrid(j,l)==1)  then ! points to be modified: sponge or land (for connectivity)
			lb = 1+modulo((j-1)-1,mgrid)
			rb = 1+modulo((j+1)-1,mgrid)
			bb = max(l-1,1)
			ub = min(l+1,ngrid)
	        datagridnew(j,l) = omega * 0.25d0 * ( datagrid(lb,l) + datagrid(rb,l) + datagrid(j,bb) + datagrid(j,ub) ) + &
	        			  (1-omega)*datagrid(j,l)
	        errij = datagridnew(j,l)-datagrid(j,l)
	        error = error + errij*errij
	    endif
      end do
    end do
    error = sqrt(error)!/count(maskgrid==1)

!    if(mod(iter,50).eq.0 ) write(*,'("cgrid=",i5,", itn=",i5,", omega=", f5.2,", err=", e10.5)'), cgrid, iter, omega, error

    do j=1,mgrid
      do l=1,ngrid
      	if (maskgrid(j,l)==1)  then
			datagrid(j,l) = datagridnew(j,l)
		endif
      end do
    end do

    iter = iter + 1
!	omega = 1/(1 - 0.25*rhojsq*omega)
  end do
!call save_matrix(datagrid, dir_grid // 'itm/' // 'c2n_hg_cgrid_after.dat')
!if (cgrid==n_grids-3) then
!	stop
!endif

  deallocate (datagrid, maskgrid)

enddo

where (mask==1) data = datagridnew

deallocate(datagridnew)

end subroutine laplace_relax


subroutine front_fill_in(data, mask)
!*****************************************************************************80
!
! FRONT_FILL_IN uses the (external!) boundary values of a connected grid-region
!	 and copies them into the interior. Implemented as a propagating boundary front

  implicit none

  integer, intent(in)		:: mask(:,:)
  real(wp), intent(inout)	:: data(:,:)

  integer :: n_tofix, n_fixed
  integer :: lb, rb, bb, ub
  integer :: j, l, dims(2)

  if (maxval( abs(shape(mask)-shape(data)) )/=0) then
  	print *, 'Error in front_fill_in:'
  	print *, 'Data and mask arrays must be of the same size'
  endif

  dims = shape(mask)



		where (mask==1) data = -1 ! mark the points that need to be modified
			n_tofix = count((data<0).and.(mask==1))
			n_fixed = 1
!			print *, n_tofix

		do while (n_fixed > 0)

			do j = 1, dims(1)
			    do l = 1, dims(2)
			    	if ((data(j,l)<0).and.(mask(j,l)==1)) then
				    	lb = max(j-1,1)
				    	rb = min(j+1,dims(1))
				    	bb = max(l-1,1)
				    	ub = min(l+1,dims(2))
				    	if (maxval(data(lb:rb,bb:ub))>0) then
				    		data(j,l) = sum( data(lb:rb,bb:ub), mask=data(lb:rb,bb:ub)>0 )/count(data(lb:rb,bb:ub)>0)
				    	endif
				    endif
				enddo
			enddo

			n_fixed = n_tofix - count((data<0).and.(mask==1)) ! how many were removed
			n_tofix = n_tofix - n_fixed
!			print *, n_tofix

		end do

		where (data<0) data = 1 ! enclosed seas/lakes filled with sponge. does not matter what value to use there

end subroutine

subroutine fd_weights(xi, x, m, dif )
!*****************************************************************************80
!
! DATA_TO_DIF calculates the pseudospectral weights c to be applied....
!
!	 Xi -- coord of the point where the approximation is calculated
!
!    Input, real(wp) X(N), the X values at which data was taken. These values must be distinct.
!
!	 M -- the order of the max derivative to be approximated. I.e. 0 -- func itself, 2 -- 2nd deriv
!
!    Output, real(wp) DIF(n,m+1), the FD coefficients for approximating the func and its M derivatives at Xi
!
  implicit none

  real(wp), intent(in)	::  x(:), xi
  real(wp), allocatable	:: dif(:,:)

  integer :: n, m, mn

  real(wp) c1, c2, c3, c4, c5
  integer i, j, k, istat


	n=size(x)-1

  if ( .not. vec_distinct( n+1, x ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
    write ( *, '(a)' ) '  Two entries of X are equal!'
    stop
  end if

!  Compute the pseudospectral weights:
!
	if (allocated(dif)) deallocate(dif)
	allocate(dif(n+1,m+1), stat = istat)
	dif = 0
	dif(1,1)=1

	c1=1
	c4=x(1)-xi

	do i = 1, n
	  mn=min(i,m)
	  c2=1
	  c5=c4
	  c4=x(i+1)-xi
	  do j = 0, i-1
	    c3=x(i+1)-x(j+1)
	    c2=c2*c3
	    if (j == i-1) then
	      do k = mn, 1, -1
	        dif(i+1,k+1)=c1*(k*dif(i,k)-c5*dif(i,k+1))/c2
	      enddo
	      dif(i+1,1)=-c1*c5*dif(i,1)/c2
	    endif

	    do k = mn, 1, -1
	      dif(j+1,k+1)=(c4*dif(j+1,k+1)-k*dif(j+1,k))/c3
	    enddo
	    dif(j+1,1)=c4*dif(j+1,1)/c3
	  enddo
	  c1=c2
	enddo

  return

end subroutine


subroutine cheby_t_zero ( n, z )
!*****************************************************************************
!
!! CHEBY_T_ZERO returns zeroes of the Chebyshev polynomial T(N)(X).
!
!  Discussion:
!
!    The I-th zero of T(N)(X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Output, real(wp) Z(N), the zeroes of T(N)(X).
!
  implicit none

  integer n

  real(wp) angle
  integer i
  real(wp), parameter :: pi = 3.141592653589793D+00
  real(wp) z(n)

  do i = 1, n
    angle = real ( 2 * i - 1,wp) * pi &
          / real ( 2 * n,    wp)
    z(i) = cos ( angle )
  end do

  return
end subroutine

subroutine cheby_u_zero ( n, z )

!*****************************************************************************
!
!! CHEBY_U_ZERO returns zeroes of the Chebyshev polynomial U(N)(X).
!
!  Discussion:
!
!    The I-th zero of U(N)(X) is cos((I-1)*PI/(N-1)), I = 1 to N
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Output, real(wp) Z(N), the zeroes of U(N)(X).
!
  implicit none

  integer n

  real(wp) angle
  integer i
  real(wp), parameter :: pi = 3.141592653589793D+00
  real(wp) z(n)

  do i = 1, n
    angle = real ( i,    wp) * pi &
          / real ( n + 1,wp)
    z(i) = cos ( angle )
  end do

  return
end subroutine
!*****************************************************************************80

function vec_distinct ( n, x )
!*****************************************************************************80
!
!! VEC_DISTINCT is true if the entries in an VEC are distinct.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real(wp) X(N), the vector to be checked.
!
!    Output, logical VEC_DISTINCT is .TRUE. if all N elements of X
!    are distinct.
!
  implicit none

  integer n

  integer i
  integer j
  logical vec_distinct
  real(wp) x(n)

  vec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( x(i) == x(j) ) then
        return
      end if
    end do
  end do

  vec_distinct = .true.

  return

end function
!*****************************************************************************

end module fd
