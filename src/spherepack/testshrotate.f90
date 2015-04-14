module testshrotate

     use precisions, only: wp, cwp
     use my_trigs
	 use control
     use save_load
	 use SHRotate
	      use spherepack_iface

      implicit none 

      contains

!**********************************************************************************************************
subroutine test_sh_rotate(dirname, filename, P)
! This function was used to perform some simple tests
     type(params) :: P

     character(len=*)	:: filename
     character(len=*) :: dirname

     integer				:: j,k
     integer, allocatable   :: topo_int(:,:)
     real(dp), allocatable	:: H_sht(:,:), topo_sht(:,:)
     real(dp), allocatable	:: lon_sht(:), lat_sht(:)
     integer			    :: nlon_sht, nlat_sht, istat
	real*8				:: angles(3)

     ! Test mars topo
	integer, parameter ::	degmax = 90, max2dgrid = 181
	integer :: lmax = 90
	real*8 ::	mcilm1(2, degmax+1, degmax+1), a(degmax+1, degmax+1), b(degmax+1, degmax+1),  header(8)
	real*8 ::	grid(2*degmax, degmax+1)
	complex*16 :: mspec((degmax+1)*(degmax+2)/2)
	character*80 ::		outsh, outgrid, infile
	integer :: m,n

	logical :: which_test(4) = (/.false.,.false.,.true.,.true./)
!	logical :: which_test(4) = (/.true.,.true.,.false.,.false./)


!%==================================================================================
!	SH Analysis -> rotate SH coeffs -> rotate SH coeffs back -> SH synthesis
!	Test 1 + Test 2: test with the 'toy' topo
!%==================================================================================

if (which_test(1) .or. which_test(2)) then
!	Get Euler angles
	print*, "Input Euler Angles"
	print*, "Alpha (deg) = "
	read(*,*) angles(1)
	print*, "Beta (deg) = "
	read(*,*) angles(2)
	print*, "Gamma (deg) = "
	read(*,*) angles(3)

angles = d2r(angles)
endif

!% Test topography
	nlon_sht = 180
	nlat_sht = 91
	    allocate(lon_sht(nlon_sht), stat = istat)
	    lon_sht = (/ ( (2*pi/nlon_sht*j), j=0,nlon_sht-1 ) /)
	    allocate(lat_sht(nlat_sht), stat = istat)
	    lat_sht = (/ ( (pi/2 - pi/(nlat_sht-1)*j), j=0,nlat_sht-1 ) /)

	allocate(H_sht(nlon_sht,nlat_sht), stat = istat)
		do j=1, nlat_sht
			do k=1, nlon_sht
				!	Test on with a rombus
				if (abs(k-nlon_sht/2) + abs(j-nlat_sht/2) < nlat_sht/2 ) then
					H_sht(k,j) = 1
					else
					H_sht(k,j) = 0
				endif
				!	TEST ON A TOY TOPO: cos(ph)*cos(th)
	!				H_sht(k,j) = real(exp(i*lon_sht(:))*cos(lat_sht(j)))
			enddo
		end do

if (which_test(1)) then
	call save_topo(trim(dirname)//'test.dat', nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht)

	call test_shrotate(H_sht, angles, (P%messages > 0), trim(dirname)//'sht_coeff_test.dat') ! spec_coeff_file
	call save_topo(trim(dirname)//'rotated_test.dat', nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht)
endif

if (which_test(2)) then
	call test_shrotate2(H_sht, angles, (P%messages > 0), trim(dirname)//'sht_coeff_test.dat') ! spec_coeff_file
	call save_topo(trim(dirname)//'rotated_test2.dat', nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht)
endif


if ( allocated(H_sht)) 	deallocate(H_sht)

!%================================================
!	SH Analysis -> SH rotate -> SH rotate back -> SH synthesis
!	Test 3 + Test 4: test with the 'real' topo
!%================================================
! Test on real topographies: SH Analysis -> rotate SH coeffs -> SH synthesis
if (which_test(3)) then
	angles =(/0, -15, -40/)
	angles = d2r(angles)

	call load_alloc_topo(trim(dirname)//'topo_rot_8.0min_pole_15_-40.dat', &
										 nlon_sht, nlat_sht, lon_sht, lat_sht, topo_int, (P%messages>=1))
	allocate(H_sht(nlon_sht,nlat_sht), stat = istat)
	H_sht = topo_int
	deallocate(topo_int)

!	print *, nlon_sht, nlat_sht, shape(H_sht)

	call test_shrotate3(H_sht, angles, (P%messages > 0),trim(dirname)//'sht_coeff_topo_rot_back_8.0min_pole_15_-40.dat') ! spec_coeff_file_out
	call save_topo(trim(dirname)//'topo_rot_back_8.0min_pole_15_-40.dat', nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht)
	deallocate(lon_sht, lat_sht, H_sht)

endif

! Test accuracy on:  SH coeffs -> rotate SH coeffs back -> SH synthesis
!
! Run which_test(3) + which_test(4)
if (which_test(4)) then
	angles =(/40, 15, 0/)
	angles = d2r(angles)

	call load_alloc_topo(trim(dirname)//'topo_rot_back_8.0min_pole_15_-40.dat', &
										 nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht, (P%messages>=1))
!	print *, nlon_sht, nlat_sht, shape(H_sht)

	call test_shrotate4(H_sht, angles, (P%messages > 0), trim(dirname)//'sht_coeff_topo_rot_back_8.0min_pole_15_-40.dat', & ! spec_coeff_file_in
												   		trim(dirname)//'sht_coeff_topo_rot_back_back_8.0min_pole_15_-40.dat') ! spec_coeff_file_out
	call save_topo(trim(dirname)//'topo_rot_back_back_8.0min_pole_15_-40.dat', nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht)
	deallocate(lon_sht, lat_sht, H_sht)

endif

end subroutine test_sh_rotate

!**********************************************************************************************************
      subroutine test_shrotate(datagrid, angles, messages, spec_coeff_file)!, datagrid_out)

      real(dp), dimension(:,:), intent(inout) :: datagrid
      character(*), intent(in), optional :: spec_coeff_file
!      real(wp), dimension(:,:), intent(inout), optional :: datagrid_out
      complex(dp), allocatable :: dataspec(:)

      real, dimension(size(datagrid,2),size(datagrid,2)) :: a, b

      integer ntrunc, n, m, nn, j, istat
      logical :: messages
	 !real(wp), allocatable	:: ph_sht(:), th_sht(:)

	real*8, intent(in)	:: angles(3)
	real*8  :: cilm1(2, size(datagrid,2),size(datagrid,2)), cilm2(2, size(datagrid,2),size(datagrid,2))
	real*8  :: dj(size(datagrid,2),size(datagrid,2),size(datagrid,2))


! perform spherical harmonic analysis.
		ntrunc = size(datagrid,2)-1
    allocate(dataspec((ntrunc+1)*(ntrunc+2)/2), stat = istat)
	call grdtospec_reg(datagrid, dataspec,'0')

      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(j,j=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         a(m,n) = 2.*real(dataspec(nn), kind=wp)
         b(m,n) = 2.*imag(dataspec(nn))
      enddo
      enddo
! analysis.

!	  call save_vector(dataspec, spec_coeff_file)
!	  deallocate(dataspec)
	  ! ------------------------------------------------------------------------------------------------

!	ROTATE SH COEFFS
	call djpi2(dj, ntrunc)	! Create rotation matrix used in the rotation routine.

	cilm1(1, :, :) = transpose(a(:,:))
	cilm1(2, :, :) = transpose(b(:,:))

	call SHRotateRealCoef(cilm2, cilm1, ntrunc, angles, dj)

	a(:, :) = transpose(cilm2(1,:,:))
	b(:, :) = transpose(cilm2(2,:,:))
      ! ------------------------------------------------------------------------------------------------

	dataspec = 0.
    dataspec = cmplx( 0.5*(/((a(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((b(m,n),n=m,ntrunc+1),m=1,ntrunc+1)/) , kind=cwp)

	call spectogrd_reg(dataspec,datagrid,'0')

end subroutine test_shrotate
!**********************************************************************************************************
subroutine test_shrotate2(datagrid, angles, messages, spec_coeff_file)!, datagrid_out)

      real(dp), dimension(:,:), intent(inout) :: datagrid
      character(*), intent(in), optional :: spec_coeff_file
!      real(wp), dimension(:,:), intent(inout), optional :: datagrid_out
      complex(dp), allocatable :: dataspec(:)

!      real, dimension(size(datagrid,2),size(datagrid,2)) :: a, b

      integer ntrunc, n, m, nn, j, istat, k
      logical :: messages
	 !real(wp), allocatable	:: ph_sht(:), th_sht(:)
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max

	real*8, intent(in)	:: angles(3)
	real*8, allocatable ::	cilm1(:,:,:), cilm2(:,:,:)
!	real*8  :: dj(size(datagrid,2),size(datagrid,2),size(datagrid,2))

! perform spherical harmonic analysis.
		ntrunc = size(datagrid,2)-1
	allocate(cilm1(2,ntrunc+1, ntrunc+1), stat = istat)
	allocate(cilm2(2,ntrunc+1, ntrunc+1), stat = istat)

    allocate(dataspec((ntrunc+1)*(ntrunc+2)/2), stat = istat)

	call grdtospec_reg(datagrid, dataspec,'0')

      cilm1 = 0.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(k,k=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         cilm1(1,n,m) = 2.*real(dataspec(nn), kind=wp)
         cilm1(2,n,m) = 2.*imag(dataspec(nn))
      enddo
      enddo

! ------------------------------------------------------------------------------------------------

!	ROTATE SH COEFFS

	call SHRotateRealCoef_2(cilm2, cilm1, ntrunc, angles)
      dataspec = cmplx( 0.5*(/((cilm2(1,n,m),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((cilm2(2,n,m),n=m,ntrunc+1),m=1,ntrunc+1)/) , kind=cwp)

!	call save_matrix_3d(cilm2, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/in/topo/nocs/cilm_new.dat')

	call spectogrd_reg(dataspec,datagrid,'0')

end subroutine test_shrotate2
!**********************************************************************************************************
subroutine test_shrotate3(datagrid, angles, messages, spec_coeff_file_out)

      real(dp), dimension(:,:), allocatable :: datagrid
      character(*), intent(in), optional :: spec_coeff_file_out
!      real(wp), dimension(:,:), intent(inout), optional :: datagrid_out
      complex(dp), allocatable :: dataspec(:)

!      real, dimension(size(datagrid,2),size(datagrid,2)) :: a, b

      integer ntrunc, n, m, nn, j, istat, k
      logical :: messages
	 !real(wp), allocatable	:: ph_sht(:), th_sht(:)
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max

	real*8, intent(in)	:: angles(3)
	real*8, allocatable ::	cilm1(:,:,:), cilm2(:,:,:)
!	real*8  :: dj(size(datagrid,2),size(datagrid,2),size(datagrid,2))

! perform spherical harmonic analysis.
	if (messages) write(*, '(a)') '1) SH analysis... '
	ntrunc = size(datagrid,2)-1
    allocate(dataspec((ntrunc+1)*(ntrunc+2)/2), stat = istat)

	call grdtospec_reg(datagrid, dataspec,'0')

	allocate(cilm1(2,ntrunc+1, ntrunc+1), stat = istat)
	allocate(cilm2(2,ntrunc+1, ntrunc+1), stat = istat)
      cilm1 = 0.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(k,k=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         cilm1(1,n,m) = 2.*real(dataspec(nn), kind=wp)
         cilm1(2,n,m) = 2.*imag(dataspec(nn))
      enddo
      enddo
! ------------------------------------------------------------------------------------------------

!	ROTATE SH COEFFS

	if (messages) write(*, '(a)') '2) rotate sht coeffs... '
	call SHRotateRealCoef_2(cilm2, cilm1, ntrunc, angles)
      dataspec = cmplx( 0.5*(/((cilm2(1,n,m),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((cilm2(2,n,m),n=m,ntrunc+1),m=1,ntrunc+1)/) , kind=cwp)
	deallocate(cilm1, cilm2)
	call save_vector(dataspec, spec_coeff_file_out)

!	call save_matrix_3d(cilm2, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/in/topo/nocs/cilm_new.dat')
	if (messages) write(*, '(a)') '3) sh synthesis... '
	call spectogrd_reg_n2(dataspec,datagrid,'0')

end subroutine test_shrotate3
!**********************************************************************************************************
subroutine test_shrotate4(datagrid, angles, messages, spec_coeff_file_in, spec_coeff_file_out)

      real(dp), dimension(:,:), allocatable :: datagrid
      character(*), intent(in), optional :: spec_coeff_file_in, spec_coeff_file_out
!      real(wp), dimension(:,:), intent(inout), optional :: datagrid_out
      complex(dp), allocatable :: dataspec(:)

!      real, dimension(size(datagrid,2),size(datagrid,2)) :: a, b

      integer ntrunc, n, m, nn, j, istat, k
      logical :: messages
	 !real(wp), allocatable	:: ph_sht(:), th_sht(:)
     real    ::      T1, T2
     integer :: wall_t1, wall_t2, clock_rate, clock_max

	real*8, intent(in)	:: angles(3)
	real*8, allocatable ::	cilm1(:,:,:), cilm2(:,:,:)
!	real*8  :: dj(size(datagrid,2),size(datagrid,2),size(datagrid,2))

! perform spherical harmonic analysis.
	  if (messages) write(*, '(a)') '1) loading SH coeffs (fast)... '
	  call load_alloc_vector(dataspec, spec_coeff_file_in)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1

	if (ntrunc /= size(datagrid,2)-1) then
		print *, "Discrepancy between the sh_coeff file and topo file dimensions:", ntrunc, size(datagrid,2)-1
		stop
	endif


	allocate(cilm1(2,ntrunc+1, ntrunc+1), stat = istat)
	allocate(cilm2(2,ntrunc+1, ntrunc+1), stat = istat)

      cilm1 = 0.
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nn = sum((/(k,k=ntrunc+1,ntrunc-m+3,-1)/))+n-m+1
         cilm1(1,n,m) = 2.*real(dataspec(nn), kind=wp)
         cilm1(2,n,m) = 2.*imag(dataspec(nn))
      enddo
      enddo
! ------------------------------------------------------------------------------------------------

!	ROTATE SH COEFFS

	if (messages) write(*, '(a)') '2) rotate sht coeffs... '
	call SHRotateRealCoef_2(cilm2, cilm1, ntrunc, angles)
      dataspec = cmplx( 0.5*(/((cilm2(1,n,m),n=m,ntrunc+1),m=1,ntrunc+1)/), &
                        0.5*(/((cilm2(2,n,m),n=m,ntrunc+1),m=1,ntrunc+1)/) , kind=cwp)
	deallocate(cilm1, cilm2)
	call save_vector(dataspec, spec_coeff_file_out)

!	call save_matrix_3d(cilm2, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/in/topo/nocs/cilm_new.dat')
	if (messages) write(*, '(a)') '3) sh synthesis... '
	call spectogrd_reg_n2(dataspec,datagrid,'0')

end subroutine test_shrotate4
!**********************************************************************************************************

! that's it!

end module testshrotate
