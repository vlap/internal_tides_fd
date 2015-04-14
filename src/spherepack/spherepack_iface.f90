module spherepack_iface

     use precisions, only: wp, cwp
     use my_trigs
     use save_load
!
! fortran95 interface for spherepack.
! For regular and gaussian lat/lon grids.
!
! requires NCAR's SPHEREPACK 3.2
! (http://www2.cisl.ucar.edu/resources/legacy/spherepack)
!
! Version 1.2 - February, 2012
! Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>

! Important Details:

! The grid and spectral arrays must be rank 2 and rank 1, respectively.
! Passing array sections is OK.

! The gridded data is assumed to be oriented such that i=1 is the 
! Greenwich meridian and j=1 is the northernmost point. Grid indices
! increase eastward and southward. If nlat is odd the equator will be
! included as a grid point. If nlat is even the equator will lie half
! way between points nlat/2 and (nlat/2)+1. nlat must be at least 3. 
! The grid increment in longitude is 2*pi/nlon radians. For example,
! nlon = 72 for a five degree grid. nlon must be greater than or 
! equal to 4. The efficiency of the computation is improved when nlon
! is a product of small prime numbers. 

! The spectral data is assumed to be in a complex array of dimension
! (MTRUNC+1)*(MTRUNC+2)/2. MTRUNC is the triangular truncation limit
! (MTRUNC = 42 for T42). MTRUNC must be <= nlat-1. Coefficients are
! ordered so that first (nm=1) is m=0,n=0, second is m=0,n=1, 
! nm=mtrunc is m=0,n=mtrunc, nm=mtrunc+1 is m=1,n=1, etc.
! In Fortran95 syntax, values of m (degree) and n (order) as a function
! of the index nm are: 

! integer, dimension((mtrunc+1)*(mtrunc+2)/2) :: indxm,indxn
! indxm = (/((m,n=m,mtrunc),m=0,mtrunc)/)
! indxn = (/((n,n=m,mtrunc),m=0,mtrunc)/)

! Conversely, the index nm as a function of m and n is: 
! nm = sum((/(i,i=mtrunc+1,mtrunc-m+2,-1)/))+n-m+1

! The associated legendre polynomials are normalized so that the
! integral (pbar(n,m,theta)**2)*sin(theta) on the interval theta=0
! to theta=pi is 1, where: 

! pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))) *
! sin(theta)**m/(2**n*factorial(n)) times the (n+m)th derivative of
! (x**2-1)**n with respect to x=cos(theta) 

! note: theta = 0.5*pi - phi, where phi is latitude and theta is colatitude.
! Therefore, cos(theta) = sin(phi) and sin(theta) = cos(phi). 
! Note that pbar(0,0,theta)=sqrt(2)/2, and 
! pbar(1,0,theta)=0.5*sqrt(6.)*sin(lat). 
!----------------------------------------------------------------------------

! The available routines are:

! SUBROUTINE GRDTOSPEC(datagrid,dataspec,gridtype):
! converts input gridded data array (datagrid) to complex spectral
! coefficients (dataspec).

! SUBROUTINE SPECTOGRD(dataspec,datagrid,gridtype):
! converts input spectral coefficient array (dataspec) to a grid (datagrid).

! SUBROUTINE SPECSMOOTH(datagrid,smooth,gridtype):
! isotropic spectral smoothing of input gridded data array (datagrid).
! Smoothing is a function of total wavenumber only, given by vector smooth
! (smooth(j),j=1,mtrunc+1, where mtrunc is the triangular truncation
! limit - see Important Details). For example, for Gaussian spectral
! smoothing, use smooth(j) = exp(-(float(j-1)/deln)**2), where deln is the
! e-folding width of the smoother in total wavenumber space. For straight
! triangular truncation, set smooth(j) = 0 for j-1 > desired truncation limit.

! SUBROUTINE CLEANUP(gridtype):
! Garbage collection. Deallocates all work arrays.
! Call this when you are done with the module to free up memory.
!----------------------------------------------------------------------------


! logical variables set to .FALSE. after work arrays initialized on
! first call. If nlon or nlat changes (as determined by size of
! input array), work arrays are real located and
! recomputed.

      implicit none 
      private 
      public :: &
      grdtospec_reg,spectogrd_reg,specsmooth_reg,cleanup_reg,&
      grdtospec_gau,spectogrd_gau,specsmooth_gau,cleanup_gau,&
      grdtospec,spectogrd,specsmooth,cleanup, &
      spectogrd_reg_n2

      integer, save :: saved_nlon_reg_spectogrd = -1  
      integer, save :: saved_nlat_reg_spectogrd = -1  
      integer, save :: saved_nlon_reg_grdtospec = -1  
      integer, save :: saved_nlat_reg_grdtospec = -1  
      integer, save :: saved_nlon_gau_spectogrd = -1  
      integer, save :: saved_nlat_gau_spectogrd = -1  
      integer, save :: saved_nlon_gau_grdtospec = -1  
      integer, save :: saved_nlat_gau_grdtospec = -1  
      logical, save :: lfrst_grdtospec_reg=.TRUE.
      logical, save :: lfrst_spectogrd_reg=.TRUE.
      logical, save :: lfrst_grdtospec_gau=.TRUE.
      logical, save :: lfrst_spectogrd_gau=.TRUE.

! work arrays are allocated in module subroutines when above
! logical variables are .TRUE., or when nlon or nlat change.

      real, dimension(:), allocatable, save :: wshaes,wshses, wshaec,wshsec!,wvhaes,wvhses
      real, dimension(:), allocatable, save :: wshags,wshsgs!,wvhags,wvhsgs

      contains
!**********************************************************************************************************
      subroutine grdtospec(datagrid,dataspec,save_sht,gridtype)
      real(wp), dimension(:,:), intent(in) :: datagrid
      complex(cwp), dimension(:), intent(out) :: dataspec
      character(*), intent(in), optional :: gridtype
      character(*), intent(in) :: save_sht

      if (present(gridtype)) then
      if (gridtype .eq. 'GAU') then
           call  grdtospec_gau(datagrid,dataspec)!save_sht not implemented
      else if (gridtype .eq. 'REG') then
           call  grdtospec_reg(datagrid,dataspec,save_sht)
      else
           write(6,*) 'SUBROUTINE GRDTOSPEC:'
           write(6,*) 'optional argument gridtype = ',gridtype
           write(6,*) 'must be either REG or GAU (default REG)'
           stop
      end if
      else
           call  grdtospec_reg(datagrid,dataspec,save_sht)
      endif
      end subroutine grdtospec
!**********************************************************************************************************
      subroutine spectogrd(dataspec,datagrid,save_sht,gridtype)
      real(wp), dimension(:,:), intent(out) :: datagrid
      complex(cwp), dimension(:), intent(in) :: dataspec
      character(*), intent(in), optional :: gridtype
      character(*), intent(in) :: save_sht

      if (present(gridtype)) then
      if (gridtype .eq. 'GAU') then
           call  spectogrd_gau(dataspec,datagrid)!save_sht not implemented
      else if (gridtype .eq. 'REG') then
            call  spectogrd_reg(dataspec,datagrid,save_sht)
      else
           write(6,*) 'SUBROUTINE SPECTOGRD:'
           write(6,*) 'optional argument gridtype = ',gridtype
           write(6,*) 'must be either REG or GAU (default REG)'
           stop
      end if
      else
           call  spectogrd_reg(dataspec,datagrid,save_sht)
      endif
      end subroutine spectogrd
!**********************************************************************************************************
      subroutine specsmooth(datagrid,smooth,gridtype, messages, spec_coeff_file)
      real(dp), dimension(:,:), intent(inout) :: datagrid
      real(dp), dimension(:), intent(in) :: smooth
      character(*), intent(in) :: gridtype
      character(*), intent(in), optional :: spec_coeff_file
      logical :: messages
!      real(wp), dimension(:,:), intent(inout), optional :: datagrid_out

      if (gridtype .eq. 'GAU') then
           call  specsmooth_gau(datagrid,smooth)!spec_coeff_file not implemented
      else if (gridtype .eq. 'REG') then
      	if (present(spec_coeff_file)) then
            call specsmooth_reg(datagrid,smooth, messages, spec_coeff_file)
	    else
	       	call specsmooth_reg(datagrid,smooth, messages)
	    endif
      else
           write(6,*) 'SUBROUTINE SPECSMOOTH:'
           write(6,*) 'optional argument gridtype = ',gridtype
           write(6,*) 'must be either REG or GAU (default REG)'
           stop
      end if
      end subroutine specsmooth
!**********************************************************************************************************
      subroutine cleanup(gridtype)
      character(*), intent(in), optional :: gridtype
      if (present(gridtype)) then
      if (gridtype .eq. 'GAU') then
           call cleanup_gau
      else if (gridtype .eq. 'REG') then
           call cleanup_reg
      else
           write(6,*) 'SUBROUTINE CLEANUP:'
           write(6,*) 'optional argument gridtype = ',gridtype
           write(6,*) 'must be either REG or GAU (default REG)'
           stop
      end if
      else
           call cleanup_reg
      endif
      end subroutine cleanup
!**********************************************************************************************************
      subroutine grdtospec_reg(datagrid,dataspec,save_sht)

! converts gridded input array (datagrid) to complex spectral coefficients
! (dataspec).

      real(wp), dimension(:,:), intent(in) :: datagrid
      complex(cwp), dimension(:), intent(out) :: dataspec

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension(2*(size(datagrid,2)+1)) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a,b
      integer nlon,nlat,lwork,ldwork,ntrunc,l1,l2,lshaes,ierror,m,n

      character(*), intent(in) :: save_sht
      logical :: load_sht_file = .false.

! compute array dimensions and infer truncation limit
! from size of spectral arrays.


      nlon = size(datagrid,1)
      nlat = size(datagrid,2)

      if (nlon .ne. saved_nlon_reg_grdtospec .or. nlat .ne. saved_nlat_reg_grdtospec) then
          lfrst_grdtospec_reg = .TRUE.
          saved_nlon_reg_grdtospec = nlon
          saved_nlat_reg_grdtospec = nlat
      end if
      lwork = size(work)
      ldwork = size(dwork)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1


! initialize work array wshaes for spherical harmonic analysis.
! only done when lfrst_grdtospec_reg = .T.
    if (lfrst_grdtospec_reg) then ! wshaes is NOT in the memory

	    ! initialize work arrays for spherical harmonic analysis.
		if (save_sht .ne. '0') then
			inquire( file=save_sht // 'sht_work.dat', exist=load_sht_file )
		else
			load_sht_file = .false.
		endif

		if ( .not. load_sht_file ) then

	      if (mod(nlon,2) .eq. 0) then
	         l1 = min0(nlat,(nlon+2)/2)
	      else
	         l1 = min0(nlat,(nlon+1)/2)
	      end if
	      if (mod(nlat,2) .eq. 0) then
	         l2 = nlat/2
	      else
	         l2 = (nlat+1)/2
	      end if

	      lshaes = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15

		  if (allocated(wshaes)) deallocate(wshaes)
	      allocate(wshaes(lshaes))
	      call shaesi(nlat,nlon,wshaes,lshaes, &
			          work,lwork,dwork,ldwork,ierror)

	      if(ierror .ne. 0) write(*,1001) ierror
 1001 	  	format(' error',i4,' in shaesi')

		!	Save to file if not kept in memory
			if (save_sht .ne. '0') then
!			    write (*, '("save wshaes ")', advance='no')
!				print *, "size(wshaes)", size(wshaes)
				call save_vector(wshaes, save_sht // 'sht_work.dat')
			endif
        else
			!	Load from file if not kept in memory
!			write (*, '("load wshaes ")', advance='no')
       	    call load_alloc_vector(wshaes, save_sht // 'sht_work.dat',lshaes)
!			lshaes = size(wshaes,1)
	    endif
	! Mark that wshaes is in memory now
	    lfrst_grdtospec_reg = .FALSE.
     else ! otherwise wshaes is in the memory
!        write (*, '("wshaes in memory ")', advance='no')
		lshaes = size(wshaes,1)
	 endif
!print *, "nlat",nlat, "nlon", nlon, "lshaes", lshaes
! transpose data.
      temp = transpose(datagrid)

! spherical harmonic analysis.
      call shaes(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshaes,lshaes,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1003) ierror
 1003  format(' error',i4,' in shaes')

! Deallocate after the loop in cpts is finished. Done in subroutine calc_sal
!		if (save_sht .ne. '0') then
!			deallocate(wshaes)
!		endif

! fill complex array dataspec with result.

      dataspec = cmplx( 0.5*(/((a(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((b(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/) , kind=cwp)
 
      end subroutine grdtospec_reg
!**********************************************************************************************************
      subroutine grdtospec_gau(datagrid,dataspec)

! converts gridded input array (datagrid) to complex spectral coefficients
! (dataspec).

      real(wp), dimension(:,:), intent(in) :: datagrid
      complex(cwp), dimension(:), intent(out) :: dataspec

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension((3*size(datagrid,2)*(size(datagrid,2)+3)+2)/2) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a,b

      integer nlon,nlat,lwork,ldwork,ntrunc,l1,l2,lshags,ierror,m,n

! compute array dimensions and infer truncation limit
! from size of spectral arrays.


      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      if (nlon .ne. saved_nlon_gau_grdtospec .or. nlat .ne. saved_nlat_gau_grdtospec) then
          lfrst_grdtospec_gau = .TRUE.
          saved_nlon_gau_grdtospec = nlon
          saved_nlat_gau_grdtospec = nlat
      end if
      lwork = size(work)
      ldwork = size(dwork)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1


! initialize work array wshags for spherical harmonic analysis.
! only done when lfrst_grdtospec_gau = .T.

      if (lfrst_grdtospec_gau) then

      if (allocated(wshags)) deallocate(wshags)

      if (mod(nlon,2) .eq. 0) then
         l1 = min0(nlat,(nlon+2)/2) 
      else
         l1 = min0(nlat,(nlon+1)/2) 
      end if
      if (mod(nlat,2) .eq. 0) then
         l2 = nlat/2      
      else
         l2 = (nlat+1)/2
      end if

      lshags = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

      allocate(wshags(lshags))

      call shagsi(nlat,nlon,wshags,lshags, &
          work,lwork,dwork,ldwork,ierror)
      if(ierror .ne. 0) write(*,1001) ierror
 1001 format(' error',i4,' in shagsi')
      lfrst_grdtospec_gau = .FALSE.
      else
      lshags = size(wshags,1)
      endif

! transpose data.

      temp = transpose(datagrid)

! spherical harmonic analysis.

      call shags(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshags,lshags,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1003) ierror
 1003  format(' error',i4,' in shags')
 
! fill complex array dataspec with result.

      dataspec = cmplx( 0.5*(/((a(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((b(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/) , kind=cwp)
 
      end subroutine grdtospec_gau
!**********************************************************************************************************
      subroutine spectogrd_reg(dataspec,datagrid,save_sht)
      
! converts complex spectral coefficients (dataspec) to 
! gridded data array (datagrid).


      real(wp), dimension(:,:), intent(inout) :: datagrid
      complex(cwp), dimension(:), intent(in) :: dataspec

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension(2*(size(datagrid,2)+1)) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a,b

      integer nlon,nlat,lwork,ldwork,ntrunc,l1,l2,lshses,ierror,m,n,nn,i

      character(*), intent(in) :: save_sht
      logical :: load_sht_file = .false.
      
! compute array dimensions and infer truncation limit
! from size of dataspec.

      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      if (nlon .ne. saved_nlon_reg_spectogrd .or. nlat .ne. saved_nlat_reg_spectogrd) then
          lfrst_spectogrd_reg = .TRUE.
          saved_nlon_reg_spectogrd = nlon
          saved_nlat_reg_spectogrd = nlat
      else
      end if
      lwork = size(work)
      ldwork = size(dwork)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1

! compute work array wshses for spherical harmonic synthesis.
! only done when lfrst_spectogrd_reg = .T.

   if (lfrst_spectogrd_reg) then ! wshses is NOT in the memory

	    ! initialize work arrays for spherical harmonic analysis.
		if (save_sht .ne. '0') then
			inquire( file=save_sht // 'shti_work.dat', exist=load_sht_file )
		else
			load_sht_file = .false.
		endif

	  if ( .not. load_sht_file ) then


	      if (mod(nlon,2) .eq. 0) then
	         l1 = min0(nlat,(nlon+2)/2)
	      else
	         l1 = min0(nlat,(nlon+1)/2)
	      end if

	      if (mod(nlat,2) .eq. 0) then
	         l2 = nlat/2
	      else
	         l2 = (nlat+1)/2
	      end if

	      lshses = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15

	      if (allocated(wshses)) deallocate(wshses)
	      allocate(wshses(lshses))

	      call shsesi(nlat,nlon,wshses,lshses,work,lwork, &
	                  dwork,ldwork,ierror)
	      if(ierror .ne. 0) write(*,1001) ierror
 1001 	format(' error',i4,' in shsesi')

		!	Save to file if not kept in memory
  		  if (save_sht .ne. '0') then
!  		    write (*, '("save wshses ")', advance='no')
			call save_vector(wshses, save_sht // 'shti_work.dat')
		  endif
       else
			!	Load from file if not kept in memory
!			write (*, '("load wshses ")', advance='no')
       	    call load_alloc_vector(wshses, save_sht // 'shti_work.dat', lshses)
!			lshses = size(wshses,1)
	    endif
	! Mark that wshses is in memory now
	    lfrst_spectogrd_reg = .FALSE.
    else ! otherwise wshses is in the memory
!        write (*, '("wshses in memory ")', advance='no')
		lshses = size(wshses,1)
	endif

! fill two real arrays (a,b) with contents of dataspec.

      a = 0.
      b = 0.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(i,i=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         a(m,n) = 2.*real(dataspec(nn), kind=wp)
         b(m,n) = 2.*imag(dataspec(nn))
      enddo
      enddo

! spherical harmonic synthesis.

      call shses(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshses,lshses,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1003) ierror
 1003  format(' error',i4,' in shses')

! Deallocate after the loop in cpts is finished. Done in subroutine calc_sal
!		if (save_sht .ne. '0') then
!			deallocate(wshses)
!		endif

! transpose data.

      datagrid = transpose(temp)
 
      end subroutine spectogrd_reg
!**********************************************************************************************************
      subroutine spectogrd_gau(dataspec,datagrid)
      
! converts complex spectral coefficients (dataspec) to 
! gridded data array (datagrid).


      real(wp), dimension(:,:), intent(out) :: datagrid
      complex(cwp), dimension(:), intent(in) :: dataspec

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension((3*size(datagrid,2)*(size(datagrid,2)+3)+2)/2) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a,b

      integer nlon,nlat,lwork,ldwork,ntrunc,l1,l2,lshsgs,ierror,m,n,nn,i
      
! compute array dimensions and infer truncation limit
! from size of dataspec.

      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      if (nlon .ne. saved_nlon_gau_spectogrd .or. nlat .ne. saved_nlat_gau_spectogrd) then
          lfrst_spectogrd_gau = .TRUE.
          saved_nlon_gau_spectogrd = nlon
          saved_nlat_gau_spectogrd = nlat
      end if
      lwork = size(work)
      ldwork = size(dwork)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1

! compute work array wshsgs for spherical harmonic synthesis.
! only done when lfrst_spectogrd_gau = .T.

      if (lfrst_spectogrd_gau) then

      if (allocated(wshsgs)) deallocate(wshsgs)

      if (mod(nlon,2) .eq. 0) then
         l1 = min0(nlat,(nlon+2)/2) 
      else
         l1 = min0(nlat,(nlon+1)/2) 
      end if
      if (mod(nlat,2) .eq. 0) then
         l2 = nlat/2        
      else
         l2 = (nlat+1)/2    
      end if

      lshsgs = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

      allocate(wshsgs(lshsgs))

      call shsgsi(nlat,nlon,wshsgs,lshsgs,work,lwork, &
                  dwork,ldwork,ierror)
      if(ierror .ne. 0) write(*,1001) ierror
 1001 format(' error',i4,' in shsgsi')
      lfrst_spectogrd_gau = .FALSE.
      else
      lshsgs = size(wshsgs,1)
      endif

! fill two real arrays (a,b) with contents of dataspec.

      a = 0.
      b = 0.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(i,i=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         a(m,n) = 2.*real(dataspec(nn), kind=wp)
         b(m,n) = 2.*imag(dataspec(nn))
      enddo
      enddo

! spherical harmonic synthesis.

      call shsgs(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshsgs,lshsgs,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1003) ierror
 1003  format(' error',i4,' in shsgs')

! transpose data.

      datagrid = transpose(temp)
 
      end subroutine spectogrd_gau
!**********************************************************************************************************
      subroutine sht_project(datagrid,ntrunc)
      	implicit none
        real(wp), dimension(:,:), intent(inout) :: datagrid
        integer, intent(in) :: ntrunc

      	print *, "blah"
      end subroutine sht_project
!**********************************************************************************************************
      subroutine specsmooth_reg(datagrid,smooth, messages, spec_coeff_file)!, datagrid_out)

! isotropic spectral smoothing of gridded data (datagrid).
! smoothing is a function of total wavenumber only, as specified
! by (smooth(n),n=1,nlat).

      real(dp), dimension(:,:), intent(inout) :: datagrid
      real(dp), dimension(:), intent(in) :: smooth
      character(*), intent(in), optional :: spec_coeff_file
!      real(wp), dimension(:,:), intent(inout), optional :: datagrid_out
      complex(dp), allocatable :: dataspec(:)

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension(2*(size(datagrid,2)+1)) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a, b

      integer nlon,nlat,lwork,ldwork,l1,l2,ierror,j, istat,m,n, nn ! lshaes,lshses
      integer lshaec,lshsec

      integer ntrunc
      logical :: messages
	  logical :: load_coeff = .false.!, use_datagrid_out = .false.
	 !real(wp), allocatable	:: ph_sht(:), th_sht(:)
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max
!----------------------------------------------------------------------
!			IMPORTANT
!Here I switched to using slower routines shaec and shsec which use O(N^2) RAM
!instead of faster shaes and shses which use O(N^3) RAM
!----------------------------------------------------------------------

! initialize the SHT coeff arrays
	a=0.
	b=0.
! compute array dimensions.

      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      lwork = size(work)
      ldwork = size(dwork)

! initialize work arrays for spherical harmonic analysis.
	if (present(spec_coeff_file)) then
		inquire( file=spec_coeff_file, exist=load_coeff )
	endif

!***************************************************
!	PERFORMANCE: START
!***************************************************
	call CPU_Time(T1)
	call system_clock ( wall_t1, clock_rate, clock_max )

	if ( .not. load_coeff ) then
		if (messages) write(*, '(a)', advance='no') '1) spherical harmonic analysis (slow)... (a file with SH coeffs is not found)'
	  ! ------------------------------------------------------------------------------------------------
      if (mod(nlon,2) .eq. 0) then
         l1 = min0(nlat,(nlon+2)/2)
      else
         l1 = min0(nlat,(nlon+1)/2)
      end if
      if (mod(nlat,2) .eq. 0) then
         l2 = nlat/2
      else
         l2 = (nlat+1)/2
      end if

!      lshaes = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
	  lshaec = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15

!      allocate(wshaes(lshaes))
      if (allocated(wshaes)) deallocate(wshaes)
      allocate(wshaec(lshaec))

!      call shaesi(nlat,nlon,wshaes,lshaes, &
!          work,lwork,dwork,ldwork,ierror)
      call shaeci(nlat,nlon,wshaec,lshaec, dwork,ldwork,ierror)
      if(ierror .ne. 0) write(*,1001) ierror
 1001 format(' error',i4,' in shaesi')
      lfrst_grdtospec_reg = .FALSE.

! transpose data.

      temp = transpose(datagrid)
!	TEST ON A TOY TOPO: cos(ph)*cos(th)
!    allocate(ph_sht(nlon), stat = istat)
!    ph_sht = (/ ( (2*pi/nlon*j), j=0,nlon-1 ) /)
!    allocate(th_sht(nlat), stat = istat)
!    th_sht = (/ ( (pi/2 - pi/(nlat-1)*j), j=0,nlat-1 ) /)
!
!	do j=1, nlat
!		temp(j, :) = real(exp(i*ph_sht(:))*cos(th_sht(j)))
!	end do

! perform spherical harmonic analysis.

!      call shaes(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
!                 wshaes,lshaes,work,lwork,ierror)
      call shaec(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshaec,lshaec,work,lwork,ierror)

      if(ierror .ne. 0) write(*,1003) ierror
 1003  format(' error',i4,' in shaes')

! fill complex array dataspec with result.
	  ntrunc = nlon/2
	  allocate(dataspec((ntrunc+1)*(ntrunc+2)/2), stat = istat)
      dataspec = cmplx( 0.5*(/((a(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((b(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/) )

!      call save_matrix(real(a,dp), spec_coeff_file//'a')
!      call save_matrix(real(b,dp), spec_coeff_file//'b')
	  call save_vector(dataspec, spec_coeff_file)
	  deallocate(dataspec)
	  ! ------------------------------------------------------------------------------------------------
	else
	  ! ------------------------------------------------------------------------------------------------
	  if (messages) write(*, '(a)', advance='no') '1) loading SH coeffs (fast)... '
	  call load_alloc_vector(dataspec, spec_coeff_file)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1
      if (ntrunc .ne. nlon/2) then
      		print *, "Error in specsmooth_reg. Inconsistency between the topo file and its SHT coeff file:", ntrunc, '/=', nlon/2
      		stop
      endif
	  ! fill two real arrays (a,b) with contents of dataspec.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(j,j=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         a(m,n) = 2.*real(dataspec(nn), kind=wp)
         b(m,n) = 2.*imag(dataspec(nn))
      enddo
      enddo
      deallocate(dataspec)
      ! ------------------------------------------------------------------------------------------------
    endif

!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (messages) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
else
	write(*, '("")')
endif

! multiply spectral coefficients by smoothing factor
! (a function of degree only).
      do j=1,nlat
         a(:,j) = smooth(j)*a(:,j)
         b(:,j) = smooth(j)*b(:,j)
      enddo

!***************************************************
!	PERFORMANCE: START
!***************************************************
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

	if (messages) write(*, '(a)', advance='no') '2) spherical harmonic synthesis... '
! initialize work arrays for spherical harmonic synthesis.
      if (mod(nlon,2) .eq. 0) then
         l1 = min0(nlat,(nlon+2)/2)
      else
         l1 = min0(nlat,(nlon+1)/2)
      end if
      if (mod(nlat,2) .eq. 0) then
         l2 = nlat/2
      else
         l2 = (nlat+1)/2
      end if

!      lshses = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
      lshsec = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15

!      if (allocated(wshses)) deallocate(wshses)
!      allocate(wshses(lshses))
      if (allocated(wshsec)) deallocate(wshsec)
      allocate(wshsec(lshsec))

!      call shsesi(nlat,nlon,wshses,lshses,work,lwork, &
!        dwork,ldwork,ierror)
      call shseci(nlat,nlon,wshsec,lshsec,dwork,ldwork,ierror)

      if(ierror .ne. 0) write(*,1002) ierror
 1002 format(' error',i4,' in shsesi')
      lfrst_spectogrd_reg = .FALSE.

! perform spherical harmonic synthesis.
!      call shses(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
!                 wshses,lshses,work,lwork,ierror)
      call shsec(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshsec,lshsec,work,lwork,ierror)

      if(ierror .ne. 0) write(*,1004) ierror
 1004  format(' error',i4,' in shses')

! transpose data.

      datagrid = transpose(temp)
!***************************************************
!	PERFORMANCE: END
!***************************************************
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

if (messages) then
	call disp ('done (CPU: ' // tostring(T2-T1) //', Wall: '&
							//tostring( real(wall_t2-wall_t1)/real(clock_rate) )//')')
else
	write(*, '("")')
endif
      end subroutine specsmooth_reg
!**********************************************************************************************************
      subroutine specsmooth_gau(datagrid,smooth)

! isotropic spectral smoothing of gridded data (datagrid).
! smoothing is a function of total wavenumber only, as specified
! by (smooth(n),n=1,nlat).

      real(dp), dimension(:,:), intent(inout) :: datagrid
      real(dp), dimension(:), intent(in) :: smooth

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension((3*size(datagrid,2)*(size(datagrid,2)+3)+2)/2) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a,b

      integer nlon,nlat,lwork,ldwork,l1,l2,lshags,lshsgs,ierror,j

! compute array dimensions.


      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      if (nlon .ne. saved_nlon_gau_spectogrd .or. &
          nlat .ne. saved_nlat_gau_spectogrd .or. &
          nlon .ne. saved_nlon_gau_grdtospec .or. &
          nlat .ne. saved_nlat_gau_grdtospec) then
          lfrst_grdtospec_gau = .TRUE.
          lfrst_spectogrd_gau = .TRUE.
          saved_nlon_gau_spectogrd = nlon
          saved_nlat_gau_spectogrd = nlat
          saved_nlon_gau_grdtospec = nlon
          saved_nlat_gau_grdtospec = nlat
      end if
      lwork = size(work)
      ldwork = size(dwork)

! if not already done, initialized work arrays for spherical
! harmonic synthesis and analysis.

      if (lfrst_grdtospec_gau) then

      if (allocated(wshags)) deallocate(wshags)

      if (mod(nlon,2) .eq. 0) then
         l1 = min0(nlat,(nlon+2)/2)
      else
         l1 = min0(nlat,(nlon+1)/2)
      end if
      if (mod(nlat,2) .eq. 0) then
         l2 = nlat/2
      else
         l2 = (nlat+1)/2
      end if

      lshags = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

      allocate(wshags(lshags))

      call shagsi(nlat,nlon,wshags,lshags, &
          work,lwork,dwork,ldwork,ierror)
      if(ierror .ne. 0) write(*,1001) ierror
 1001 format(' error',i4,' in shagsi')
      lfrst_grdtospec_gau = .FALSE.
      else
      lshags = size(wshags,1)
      end if
      if (lfrst_spectogrd_gau) then

      if (allocated(wshsgs)) deallocate(wshsgs)

      if (mod(nlon,2) .eq. 0) then
         l1 = min0(nlat,(nlon+2)/2)
      else
         l1 = min0(nlat,(nlon+1)/2)
      end if
      if (mod(nlat,2) .eq. 0) then
         l2 = nlat/2
      else
         l2 = (nlat+1)/2
      end if

      lshsgs = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

      allocate(wshsgs(lshsgs))

      call shsgsi(nlat,nlon,wshsgs,lshsgs,work,lwork, &
        dwork,ldwork,ierror)
      if(ierror .ne. 0) write(*,1002) ierror
 1002 format(' error',i4,' in shsgsi')
      lfrst_spectogrd_gau = .FALSE.
      else
      lshsgs = size(wshsgs,1)
      endif

! transpose data.

      temp = transpose(datagrid)

! perform spherical harmonic analysis.

      call shags(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshags,lshags,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1003) ierror
 1003  format(' error',i4,' in shags')

! multiply spectral coefficients by smoothing factor
! (a function of degree only).

      do j=1,nlat
         a(:,j) = smooth(j)*a(:,j)
         b(:,j) = smooth(j)*b(:,j)
      enddo

! perform spherical harmonic synthesis.

      call shsgs(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshsgs,lshsgs,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1004) ierror
 1004  format(' error',i4,' in shsgs')

! transpose data.

      datagrid = transpose(temp)

      end subroutine specsmooth_gau
!**********************************************************************************************************
      subroutine cleanup_gau

! cleanup memory allocations when done with module.

      if (allocated(wshags)) deallocate(wshags)
      if (allocated(wshsgs)) deallocate(wshsgs)
!      if (allocated(wvhags)) deallocate(wvhags)
!      if (allocated(wvhsgs)) deallocate(wvhsgs)
      saved_nlat_gau_grdtospec = -1
      saved_nlon_gau_grdtospec = -1
      saved_nlat_gau_spectogrd = -1
      saved_nlon_gau_spectogrd = -1
      lfrst_grdtospec_gau = .TRUE.
      lfrst_spectogrd_gau = .TRUE.

      end subroutine cleanup_gau
!**********************************************************************************************************
      subroutine cleanup_reg

! cleanup memory allocations when done with module.

      if (allocated(wshaes)) deallocate(wshaes)
      if (allocated(wshses)) deallocate(wshses)
      if (allocated(wshaec)) deallocate(wshaec)
      if (allocated(wshsec)) deallocate(wshsec)
!      if (allocated(wvhaes)) deallocate(wvhaes)
!      if (allocated(wvhses)) deallocate(wvhses)
      saved_nlat_reg_grdtospec = -1
      saved_nlon_reg_grdtospec = -1
      saved_nlat_reg_spectogrd = -1
      saved_nlon_reg_spectogrd = -1
      lfrst_grdtospec_reg = .TRUE.
      lfrst_spectogrd_reg = .TRUE.

      end subroutine cleanup_reg

!**********************************************************************************************************
!	SPHERICAL HARMONIC ANALYSIS and SYNTHESIS routines that are using O(N**2) storage instead of O(N**3)
!**********************************************************************************************************
!----------------------------------------------------------------------
!						IMPORTANT
!	Here I switched to using slower routines shaec and shsec which use O(N^2) RAM
!	instead of faster shaes and shses which use O(N^3) RAM
!	NOT TO BE USED FOR REPETITIVE COMPUTATIONS
!----------------------------------------------------------------------

      subroutine spectogrd_reg_n2(dataspec,datagrid,save_sht)

! converts complex spectral coefficients (dataspec) to
! gridded data array (datagrid).


      real(wp), dimension(:,:), intent(inout) :: datagrid
      complex(cwp), dimension(:), intent(in) :: dataspec

      real, dimension((4*size(datagrid,1)+2)*size(datagrid,2)) :: work
      double precision, dimension(2*(size(datagrid,2)+1)) :: dwork
      real, dimension(size(datagrid,2),size(datagrid,1)) :: temp
      real, dimension(size(datagrid,2),size(datagrid,2)) :: a,b

      integer nlon,nlat,lwork,ldwork,ntrunc,l1,l2,ierror,m,n,nn,i
      integer lshsec

      character(*), intent(in) :: save_sht
      logical :: load_sht_file = .false.

! compute array dimensions and infer truncation limit
! from size of dataspec.

      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      if (nlon .ne. saved_nlon_reg_spectogrd .or. nlat .ne. saved_nlat_reg_spectogrd) then
          lfrst_spectogrd_reg = .TRUE.
          saved_nlon_reg_spectogrd = nlon
          saved_nlat_reg_spectogrd = nlat
      else
      end if
      lwork = size(work)
      ldwork = size(dwork)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1

! compute work array wshsec for spherical harmonic synthesis.
! only done when lfrst_spectogrd_reg = .T.

   if (lfrst_spectogrd_reg) then ! wshsec is NOT in the memory

	    ! initialize work arrays for spherical harmonic analysis.
		if (save_sht .ne. '0') then
			inquire( file=save_sht // 'shti_n2_work.dat', exist=load_sht_file )
		else
			load_sht_file = .false.
		endif

	  if ( .not. load_sht_file ) then


	      if (mod(nlon,2) .eq. 0) then
	         l1 = min0(nlat,(nlon+2)/2)
	      else
	         l1 = min0(nlat,(nlon+1)/2)
	      end if

	      if (mod(nlat,2) .eq. 0) then
	         l2 = nlat/2
	      else
	         l2 = (nlat+1)/2
	      end if

	      lshsec = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15

	      if (allocated(wshsec)) deallocate(wshsec)
	      allocate(wshsec(lshsec))

		  call shseci(nlat,nlon,wshsec,lshsec,dwork,ldwork,ierror)

      	  if(ierror .ne. 0) write(*,1003) ierror
 1003 	format(' error',i4,' in shseci')
 		  lfrst_spectogrd_reg = .FALSE.

		!	Save to file if not kept in memory
  		  if (save_sht .ne. '0') then
!  		    write (*, '("save wshsec ")', advance='no')
			call save_vector(wshsec, save_sht // 'shti_work.dat')
		  endif
       else
			!	Load from file if not kept in memory
!			write (*, '("load wshsec ")', advance='no')
       	    call load_alloc_vector(wshsec, save_sht // 'shti_n2_work.dat', lshsec)
!			lshsec = size(wshsec,1)
	    endif
	! Mark that wshsec is in memory now
	    lfrst_spectogrd_reg = .FALSE.
    else ! otherwise wshsec is in the memory
!        write (*, '("wshsec in memory ")', advance='no')
		lshsec = size(wshsec,1)
	endif

! fill two real arrays (a,b) with contents of dataspec.

      a = 0.
      b = 0.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(i,i=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         a(m,n) = 2.*real(dataspec(nn), kind=wp)
         b(m,n) = 2.*imag(dataspec(nn))
      enddo
      enddo

! spherical harmonic synthesis.

      call shsec(nlat,nlon,0,1,temp,nlat,nlon,a,b,nlat,nlat, &
                 wshsec,lshsec,work,lwork,ierror)
      if(ierror .ne. 0) write(*,1004) ierror
 1004  format(' error',i4,' in shsec')

! Deallocate after the loop in cpts is finished. Done in subroutine calc_sal
!		if (save_sht .ne. '0') then
!			deallocate(wshsec)
!		endif

! transpose data.

      datagrid = transpose(temp)

      end subroutine spectogrd_reg_n2
!**********************************************************************************************************

! that's it!

      end module spherepack_iface
