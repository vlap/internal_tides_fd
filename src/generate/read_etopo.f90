module read_etopo

!     use save_load
     use my_trigs
     use err_manager
     use precisions, only: dp

     implicit none

     integer, parameter               :: ilon = 21601
     integer, parameter               :: jlat = 10801

!	NOCS variables
!     character(len=*), parameter               :: x_name = 'lon'
!     character(len=*), parameter               :: y_name = 'lat'
!     character(len=*), parameter               :: z_name = 'Bathy_meters'
!	ETOPO variables
     character(len=*), parameter               :: x_name = 'x'
     character(len=*), parameter               :: y_name = 'y'
     character(len=*), parameter               :: z_name = 'z'

     real(dp), parameter  ::      eps_tmp = 1e-1                  ! epsilon for rounding error


        contains

!********************************************************************
!********************************************************************
!     MAIN SUBROUTINE TO READ IN GRID
!      (all data is in degrees).
!********************************************************************
!********************************************************************

subroutine prepare_topo(etopo_file, numLonsR, numLatsR, xValuesR, yValuesR, zValuesR, lonP, latP)
     integer::numLons, numLats
     real(dp), allocatable :: xValues(:), yValues(:)
     integer, allocatable :: zValues(:, :)

     !integer     ::     numLonsI, numLatsI
     !real(dp), allocatable :: xValuesI(:), yValuesI(:)
     !integer, allocatable :: zValuesI(:, :)

     real(dp)     ::     latP, lonP, lonS
     integer    ::     numLonsR, numLatsR
     real(dp), allocatable :: xValuesR(:), yValuesR(:)
     integer, allocatable :: zValuesR(:, :)

     character(len=*), intent(in) :: etopo_file


      call read_etopo_topo(numLons, numLats, xValues, yValues, zValues, etopo_file)
      print *, yValues(1), yValues(numLats)

     numLonsR = numLonsR + 1
     numLatsR = numLatsR + 1
     if ( (numLonsR > numLons) .or. (numLatsR > numLats) ) then
		numLonsR = numLons
		numLatsR = numLats
	 endif

!  call load_alloc_topo('/scratch/1/users/amtvl/Tides/data/NOCS/topo_2.0min.dat', numLons, numLats, xValues, yValues, zValues)
!     numLonsR = numLons;
!     numLatsR = numLats;

     call rotate_etopo(numLons, numLats, xValues, yValues, zValues, &
        latP, lonP, numLonsR, numLatsR, xValuesR, yValuesR, zValuesR)

     deallocate(xValues)
     deallocate(yValues)
     deallocate(zValues)

!      careful, lonR must match the number of grid points !!!!! lonS = 180
	lonS = 180
      call shift_etopo(numLonsR, numLatsR, xValuesR, zValuesR, lonS)

end subroutine prepare_topo

!********************************************************************
!********************************************************************
!     READ NETCDF
!********************************************************************
!********************************************************************
     subroutine read_etopo_topo(numLons, numLats, xValues, yValues, zValues, filename)

!     include 'netcdf.inc'
     use netcdf

        integer ::      ncid, status, statusx, statusy, statusz, &
                       xdimid, ydimid,                     &     ! Dimensions ID
                       xvarid, yvarid, zvarid,           &     ! Variable ID
                numDims, numAtts,                &
                actualRangeLength, titleLength, long_nameLength, unitsLength ! Attribute lengths

    real, dimension(:), allocatable               :: XactualRange, YactualRange
    character (len = 80)     ::     Xlong_name, Ylong_name
    character (len = 80)     ::     Xunits, Yunits
    character (len = 80)     ::     title
    character (len = *)      :: filename

    integer, dimension(nf90_max_var_dims) :: DimIds

     integer     ::     numLons, numLats
     real(dp), allocatable :: xValues(:), yValues(:)
     integer, allocatable :: zValues(:, :)

     write(*,'("Loading ETOPO netcdf file: ", a)', advance='no'), filename
     write(*, '(" at resolution ", f4.2, " minute: ")', advance='no')   360*60/real(ilon-1)

!     Open netcdf file
    status = nf90_open(filename, nf90_NoWrite, ncid)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Open netcdf file');
          stop
     end if

!     Inquire Dimensions ID
    status = nf90_inq_dimid(ncid, x_name, xdimid)
    if (status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Dimensions ID');
          stop
     end if

    status = nf90_inq_dimid(ncid, y_name, ydimid)
    if (status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Dimensions ID');
          stop
     end if

!     Inquire Variables ID
    status = nf90_inq_varid(ncid, x_name, xvarid)
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variables ID');
          stop
     end if

    status = nf90_inq_varid(ncid, y_name, yvarid)
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variables ID');
          stop
     end if
    status = nf90_inq_varid(ncid, z_name, zvarid)
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variables ID');
          stop
     end if

!     Inquire Variables
     !     X                                                                                                    !optional begins
      !     Inquire Attributes
    status = nf90_inquire_attribute(ncid, xvarid, "actual_range", len = actualRangeLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes X');
          stop
     end if
    status = nf90_inquire_attribute(ncid, xvarid, "long_name", len = long_nameLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes X');
          stop
     end if
     status = nf90_inquire_attribute(ncid, xvarid, "units", len = unitsLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes X');
          stop
     end if
    !Allocate space to hold attribute values
    allocate(XactualRange(actualRangeLength), stat = statusx)
    if(statusx /= 0 .or. len(Xunits) < unitsLength .or.  len(Xlong_name) < long_nameLength ) then
         print *, "Not enough space to put values of X attributes."
        stop
         end if
    ! Read the attributes.
    status = nf90_get_att(ncid, xvarid, "long_name", Xlong_name)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes X');
          stop
     end if
    status = nf90_get_att(ncid, xvarid, "actual_range", XactualRange)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes X');
          stop
     end if
    status = nf90_get_att(ncid, xvarid, "units", Xunits)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes X');
          stop
     end if
       !     Y
      !     Inquire Attributes
    status = nf90_inquire_attribute(ncid, yvarid, "actual_range", len = actualRangeLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes Y');
          stop
     end if
    status = nf90_inquire_attribute(ncid, yvarid, "long_name", len = long_nameLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes Y');
          stop
     end if
    status = nf90_inquire_attribute(ncid, yvarid, "units", len = unitsLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes Y');
          stop
     end if
    !Allocate space to hold attribute values
    allocate(YactualRange(actualRangeLength), stat = statusy)
    if(statusy /= 0 .or. len(Yunits) < unitsLength .or.  len(Ylong_name) < long_nameLength ) then
            print *, "Not enough space to put values of Y attributes."
        stop
         end if
    ! Read the attributes.
    status = nf90_get_att(ncid, Yvarid, "long_name", Ylong_name)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes Y');
          stop
     end if

    status = nf90_get_att(ncid, Yvarid, "actual_range", YactualRange)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes Y');
          stop
     end if
    status = nf90_get_att(ncid, Yvarid, "units", Yunits)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes Y');
          stop
     end if

     !     Global
    status = nf90_inquire_attribute(ncid, nf90_global, "title", len = titleLength)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Inquire Attributes');
          stop
     end if
     !Check string lengths
    if(status /= nf90_NoErr .or. len(title) < titleLength) then
         print *, "Not enough space to put attribute values."
        stop
         end if
    status = nf90_get_att(ncid, nf90_global, "title", title)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Read Attributes X');
          stop
     end if                                                                      !optional ends

     !     Z
    status = nf90_inquire_variable(ncid, zvarid, ndims = numDims, natts = numAtts)
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variable Z');
          stop
     end if
    status = nf90_inquire_variable(ncid, zvarid, dimids = DimIds(:numDims))
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variable Z');
          stop
     end if

    status = nf90_inquire_dimension(ncid, dimIDs(1), len = numLons)
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variable Z');
          stop
     end if
    status = nf90_inquire_dimension(ncid, dimIDs(2), len = numLats)
    if(status /= nf90_NoErr)  then
          call handle_nc_err(status, 'Inquire Variable Z');
          stop
     end if

!     Read Values
    if (allocated(ZValues)) deallocate(ZValues)
    allocate(ZValues(numLons, numLats), stat = statusz)
    if (allocated(XValues)) deallocate(XValues)
    allocate(XValues(numLons), stat = statusx)
    if (allocated(YValues)) deallocate(YValues)
    allocate(YValues(numLats), stat = statusy)

    if (statusx /= 0 .or. statusy /= 0 .or. statusz /= 0 .or. numLons /= ilon .or. numLats /= jlat) then
         print *, "[Topo] Can't allocate the memory."
        stop
    end if

    statusx = nf90_get_var(ncid, xVarId, xValues)
    statusy = nf90_get_var(ncid, yVarId, yValues)
    statusz = nf90_get_var(ncid, zVarId, zValues)
    if(statusx /= nf90_NoErr .or. statusy /= nf90_NoErr .or. statusz /= nf90_NoErr ) then
          call handle_nc_err(statusx, 'Reading netcdf file X');
          call handle_nc_err(statusy, 'Reading netcdf file Y');
          call handle_nc_err(statusz, 'Reading netcdf file Z');
          stop
    end if

!     Close netcdf file
    status = nf90_close(ncid)
    if (status /= nf90_noerr)  then
          call handle_nc_err(status, 'Closing netcdf file');
          stop
     end if

!     print *, "status =",status, ", ncid =", ncid



!     write(*,*) "size Z =",shape(zValues), ",size X =",size(xValues), ",size Y =",size(yValues)
!     Explore the Contents

     ![numdims,nvars,natts] = netcdf.inq(ncid);


!     call opngks                                ! open graphics device
     !call cpcnrc(zValues, numLons, numLons, numLats, 0.,0.,0.,0,0,0)  ! make the plot
     !call frame                                 ! close out the plot
     !call clsgks                                ! close graphics device


!        print *, "hi from module1"
!        print *, ilon, jlat
        print *, "done."
        print *, ""

        end subroutine
!********************************************************************
!********************************************************************
! Create the topography file, which should include
! the fields:
!  zValuesR, yValuesR, xValuesR, latp, lonp
! where topo is a lat x lon grid of height (m) above sea level,
! lat and lon are the corresponding arbitrarily spaced vectors,
! and latp and lonp give the position of the rotated pole
! (all in degrees).
!********************************************************************
!********************************************************************
    subroutine rotate_etopo(numLons, numLats, xValues, yValues, zValues, &
        latP, lonP, numLonsR, numLatsR, xValuesR, yValuesR, zValuesR)

     integer ::      statusx, statusy, statusz, statusx1, statusy1
     integer ::      k, j, indi, indj !, tmp
     real(dp)  ::      dx, dy, dxR, dyR
     real(dp), allocatable :: xMap_gb(:, :), yMap_gb(:, :)     !map the new rotated grid back to the original coordinates

     real(dp)     ::     weight1, weight2
     integer, intent(in)     ::     numLons, numLats
     real(dp), allocatable, intent(in) :: xValues(:), yValues(:)
     integer, allocatable, intent(in) :: zValues(:, :)

     real(dp)     ::     latP, lonP, ztmp1, ztmp2
     integer    ::     numLonsR, numLatsR
     real(dp), allocatable :: xValuesR(:), yValuesR(:)
     integer, allocatable :: zValuesR(:, :)

     write(*,'(a, f6.2, ", ", f6.2, ") degrees")') "Transform to lat-lon with the poles at (", lonP, 90-latP
     write(*, '("and convert to the ", f5.2, " min resolution: ")', advance='no')   360*60/real(numLonsR-1)

If (numLons < numLonsR .or. numLats < numLatsR .or. numLonsR <= 0 .or. numLatsR <= 0 ) then
     Write (*, *) 'Bad Dimenisons...subroutine rotate_etopo'
     stop
end if

If (xValues(1) >= xValues(numLons) .or. yValues(1) >= yValues(numLats) ) then
     Write (*, *) 'Bad Inputs...subroutine rotate_etopo'
     stop
end if

     if (allocated(xMap_gb)) deallocate(xMap_gb)
    allocate(xMap_gb(numLonsR, numLatsR), stat = statusx1)
     if (allocated(yMap_gb)) deallocate(yMap_gb)
    allocate(yMap_gb(numLonsR, numLatsR), stat = statusy1)
     if (allocated(ZValuesR)) deallocate(ZValuesR)
    allocate(ZValuesR(numLonsR, numLatsR), stat = statusz)
    if (allocated(XValuesR)) deallocate(XValuesR)
    allocate(XValuesR(numLonsR), stat = statusx)
    if (allocated(YValuesR)) deallocate(YValuesR)
    allocate(YValuesR(numLatsR), stat = statusy)

    if (statusx1 /= 0 .or. statusy1 /= 0 .or. statusx /= 0 .or. statusy /= 0 .or. statusz /= 0) then
         print *, "[Topo] Can't allocate the memory."
        stop
    end if

    dx = (xValues(numLons) - xValues(1))/(numLons - 1) !the topography is periodical, i.e. map from -180 to 180
    dy = (yValues(numLats) - yValues(1))/(numLats - 1) !and data for xValues(1) and xValues(numLons) is the same
    dxR = (xValues(numLons) - xValues(1))/(numLonsR - 1)
    dyR = (yValues(numLats) - yValues(1))/(numLatsR - 1)


     XValuesR = (/ ( (xValues(1) + dxR*(j-1)), j=1,numLonsR ) /)
     yValuesR = (/ ( (yValues(1) + dyR*(j-1)), j=1,numLatsR ) /)

!********************************************************************
!     treat poles separately!
!     New south pole
!print *, "Hi south pole"
         xMap_gb(1, 1) = -lonP
         yMap_gb(1, 1) = -(90 - latP)
         indi = (xMap_gb(1, 1) - xValues(1))/dx  !    convert to integer
         indj = (yMap_gb(1, 1) - yValues(1))/dy  !    convert to integer

          weight1 = (xMap_gb(1, 1) - xValues(indi+1))/dx
          weight2 = (xValues(indi+2) - xMap_gb(1, 1))/dx

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 1', weight1, weight2);
                stop
          end if

         !     First horizontal averaging
          ztmp1 = ZValues(indi + 2, indj+1)*weight1 + &
                              ZValues(indi+1, indj+1)*weight2
          ztmp2 = ZValues(indi + 2, indj+2)*weight1 + &
                              ZValues(indi+1, indj+2)*weight2

         !     Then vertical
          weight1 = (yMap_gb(1, 1) - yValues(indj+1))/dy
          weight2 = (yValues(indj+2) - yMap_gb(1, 1))/dy

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 2', weight1, weight2);
!                write(*,*) (indi(j), j=1,k+1)
!                write(*,*) ((xValuesI(j) - xValues(indi(j)))/dx, j=2,k+1)
!                print *, indj, yValues(indj), yMap_gb(1, 1), indj*dy + yValues(1), yValues(indj)
!                write(*,'(f12.10)'), 1. + eps_tmp
                stop
          end if


          ZValuesR(1, 1) = ztmp2*weight1 + ztmp1*weight2

!print *, ZValuesR(1, 1)

    do k = 2, numLonsR
         xMap_gb(k, 1) = xMap_gb(k-1, 1)
         yMap_gb(k, 1) = yMap_gb(k-1, 1)
         ZValuesR(k, 1) = ZValuesR(k-1, 1)
    end do

!********************************************************************
!     Then new north pole
!print *, "Hi north pole"
         xMap_gb(1, numLatsR) = lonP
         yMap_gb(1, numLatsR) = (90 - latP)

         indi = (xMap_gb(1, numLatsR) - xValues(1))/dx  !    convert to integer
         indj = (yMap_gb(1, numLatsR) - yValues(1))/dy  !    convert to integer

         !     First horizontal averaging
          weight1 = (xMap_gb(1, numLatsR)  - xValues(indi+1))/dx
          weight2 = (xValues(indi+2) - xMap_gb(1, numLatsR) )/dx

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 3', weight1, weight2);
                stop
          end if

          ztmp1 = ZValues(indi + 2, indj+1)*weight1 + &
                              ZValues(indi+1, indj+1)*weight2
          ztmp2 = ZValues(indi + 2, indj+2)*weight1 + &
                              ZValues(indi+1, indj+2)*weight2

         !     Then vertical
          weight1 = (yMap_gb(1, numLatsR) - yValues(indj+1))/dy
          weight2 = (yValues(indj+2) - yMap_gb(1, numLatsR))/dy

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 4', weight1, weight2);
                stop
          end if

          ZValuesR(1, numLatsR) = ztmp2*weight1 + ztmp1*weight2

!print *, ZValuesR(1, numLatsR)
    do k = 2, numLonsR
         xMap_gb(k, numLatsR) = xMap_gb(k-1, numLatsR)
         yMap_gb(k, numLatsR) = yMap_gb(k-1, numLatsR)
         ZValuesR(k, numLatsR) = ZValuesR(k-1, numLatsR)
    end do

!********************************************************************
!   And then for the bulk of the map
!print *, "Hi all"
!     First convert everything in degrees to radians
    latP = d2r(latP)
    lonP = d2r(lonP)
!    do k = 1, numLonsR
!         xValuesR(k) = d2r(xValuesR(k))
!    end do
!    do j = 1, numLatsR
!         yValuesR(j) = d2r(yValuesR(j))
!    end do
	xValuesR = d2r(xValuesR)
	yValuesR = d2r(yValuesR)

!print *, "type in something"
!read (*,*) tmp
    do k = 1, numLonsR
         do j = 2, numLatsR-1

         yMap_gb(k, j) =  asin( - sin(latP)*cos(yValuesR(j))*cos(xValuesR(k)) + &
                                   cos(latP)*sin(yValuesR(j)) )

!         xMap_gb(k, j) = lonP + asin( cos(yValuesR(j))*sin(xValuesR(k))/cos(yMap_gb(k, j)) )
         if (cos(yValuesR(j))*sin(xValuesR(k))/cos(yMap_gb(k, j)) > 0) then
              xMap_gb(k, j) = lonP + acos( max(min( ( cos(latP)*cos(yValuesR(j))*cos(xValuesR(k)) + &
                                   sin(latP)*sin(yValuesR(j)) )/cos(yMap_gb(k, j)), 1.), -1.) )
!         write (*,'("if 1, ", "k:", i4, ", j:", i4, ", x:", f10.3, ", y:", f10.3)') k, j, xMap_gb(k, j), yMap_gb(k, j)
         else
              xMap_gb(k, j) = lonP - acos( max(min( ( cos(latP)*cos(yValuesR(j))*cos(xValuesR(k)) + &
                                   sin(latP)*sin(yValuesR(j)) )/cos(yMap_gb(k, j)), 1.), -1.) )
!         write (*,'("if 2, ", "k:", i4, ", j:", i4, ", x:", f10.3, ", y:", f10.3)') k, j, xMap_gb(k, j), yMap_gb(k, j)
         end if
         if (xMap_gb(k, j) < -pi) then
              xMap_gb(k, j) = xMap_gb(k, j) + 2*pi
         elseif (xMap_gb(k, j) > pi) then
            xMap_gb(k, j) = xMap_gb(k, j) + 2*pi
         endif

         xMap_gb(k, j) =  r2d(xMap_gb(k, j))
         yMap_gb(k, j) =  r2d(yMap_gb(k, j))

!if ((k == (numLonsR-1)/2+1).and.(j == (numLatsR-1)/2+1)) then
!print *, "1", r2d(xValuesR(k)), r2d(yValuesR(j))
!print *, "2", xMap_gb(k, j)
!print *, "3", yMap_gb(k, j)
!print *, "4", min((cos(latP)*cos(yValuesR(j))*cos(xValuesR(k)) + sin(latP)*sin(yValuesR(j)) )/cos(d2r(yMap_gb(k, j))),1.)
!print *, "5",  acos( cos(latP)*cos(yValuesR(j))*cos(xValuesR(k)) + sin(latP)*sin(yValuesR(j)) )/cos(d2r(yMap_gb(k, j)))

!endif
         enddo
    enddo

!     Convert everything back to degrees
    latP = r2d(latP)
    lonP = r2d(lonP)
    xValuesR = r2d(xValuesR)
    yValuesR = r2d(yValuesR)

    do k = 1, numLonsR
         do j = 2, numLatsR-1

         indi = (xMap_gb(k, j) - xValues(1))/dx  !    convert to integer
         indj = (yMap_gb(k, j) - yValues(1))/dy  !    convert to integer

         !     First horizontal averaging
          weight1 = (xMap_gb(k, j) - xValues(indi+1))/dx
          weight2 = (xValues(indi+2) - xMap_gb(k, j))/dx

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 5', weight1, weight2);
                stop
          end if
!if ((weight1<0).or.(weight2<0)) then
!       write (*,'("k:", i4, ", j:", i4, ", indi:", i5, ", indj:", i5, ", weightx left:", f10.3, ", weightx right:", f10.3)') k, j, &
!          indi, indj, weight1, weight2
!end if
          ztmp1 = ZValues(indi + 2, indj+1)*weight1 + &
                              ZValues(indi+1, indj+1)*weight2
          ztmp2 = ZValues(indi + 2, indj+2)*weight1 + &
                              ZValues(indi+1, indj+2)*weight2
         !     Then vertical
          weight1 = (yMap_gb(k, j) - yValues(indj+1))/dy
          weight2 = (yValues(indj+2) - yMap_gb(k, j))/dy

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 6', weight1, weight2);
                stop
          end if

          ZValuesR(k, j) = ztmp2*weight1 + ztmp1*weight2

!if (abs(ZValuesR(k,j))> 10000) then
!       write (*,'("k:", i4, ", j:", i4, ", ZValuesR(k, j):", i10)') k, j, ZValuesR(k, j)
!end if
!       write (*,'("k:", i4, ", j:", i4, ", indi:", i5, ", indj:", i5, ", weighty left:", f10.3, ", weighty right:", f10.3)') k, j, &
          !indi, indj, weight1, weight2
         enddo
    enddo

!print *, MINVAL(ZValuesR), MINLOC(ZValuesR)
!print *, MAXVAL(ZValuesR), MAXLOC(ZValuesR)


     if (allocated(xMap_gb)) deallocate(xMap_gb)
     if (allocated(yMap_gb)) deallocate(yMap_gb)

        print *, "done."
        print *, ""

    end subroutine

!********************************************************************
!********************************************************************
! Shift the lat-lon topography (rotation around the pole)
! lonR gives the position rotation angle (all in degrees).
!********************************************************************
!********************************************************************
    subroutine shift_etopo(numLons, numLats, xValues, zValues, lonR)

     integer ::      shift
     integer, intent(in)     ::     numLons, numLats
     real(dp) :: xValues(:)
     integer,  intent(inout) :: zValues(:, :)
     real(dp), intent(in)      ::     lonR

     integer ::               status
     real(dp), allocatable ::       ztmp(:, :)
     !integer, allocatable ::      permut(:)

    if (allocated(ztmp)) deallocate(ztmp)
    allocate(ztmp(numLons, numLats), stat = status)
    ztmp = zValues

    !if (allocated(permut)) deallocate(permut)
!    allocate(permut(numLons), stat = status)

    xValues = xValues + lonR
!	must be integer
    shift = (numLons-1)*lonR/360
!    print *, numLons,shift

!   assumes that zValues(:, 1) = zValues(:, numLons)
    zValues(1:(numLons-shift), :) = ztmp(shift+1:numLons, :)
    zValues((numLons-shift + 1):numLons, :) = ztmp(2:shift+1, :)

    !permut = (/ ( mod(k + shift, numLons), k=1,numLons ) /)
!    do k = 1, numLons
!         zValues(:, k) = ztmp(:, permut(k))
!    enddo

    deallocate(ztmp, stat = status)

    end subroutine

!********************************************************************
!********************************************************************
!    Interpolate equal spaced grid lat-lon to lower resolution (using simple bilinear interpolation)
!********************************************************************
!********************************************************************
    subroutine interp_etopo(numLons, numLats, xValues, yValues, zValues, &
        numLonsI, numLatsI, xValuesI, yValuesI, zValuesI)


     integer ::      statusx, statusy, statusz, status
     integer ::      k, j
     integer, allocatable :: indi(:), indj(:)
     real(dp)  ::      dx, dy, dxI, dyI

     integer     ::     numLons, numLats
     real(dp), allocatable :: xValues(:), yValues(:)
     integer, allocatable :: zValues(:, :)

     real(dp)     ::     weight1, weight2
     real(dp)     ::     ztmp1, ztmp2
     integer     ::     numLonsI, numLatsI
     real(dp), allocatable :: xValuesI(:), yValuesI(:)
     integer, allocatable :: zValuesI(:, :)

If (numLons <= numLonsI .or. numLats <= numLatsI .or. numLonsI <= 0 .or. numLatsI <= 0 ) then
     Write (*, *) 'Bad Dimenisons...subroutine interp_etopo'
     stop
end if

If (xValues(1) >= xValues(numLons) .or. yValues(1) >= yValues(numLats) ) then
     Write (*, *) 'Bad Inputs...subroutine interp_etopo'
     stop
end if

	if (allocated(ZValuesI)) deallocate(ZValuesI)
    allocate(ZValuesI(numLonsI, numLatsI), stat = statusz)
    if (allocated(XValuesI)) deallocate(XValuesI)
    allocate(XValuesI(numLonsI), stat = statusx)
    if (allocated(YValuesI)) deallocate(YValuesI)
    allocate(YValuesI(numLatsI), stat = statusy)

    if (statusx /= 0 .or. statusy /= 0 .or. statusz /= 0) then
         print *, "[Topo] Can't allocate the memory."
        stop
    end if

    dx = (xValues(numLons) - xValues(1))/(numLons - 1) !the topography is periodical, k.e. map from-180 to 180
    dy = (yValues(numLats) - yValues(1))/(numLats - 1) !and data for xValues(1) and xValues(numLons) is the same
    dxI = (xValues(numLons) - xValues(1))/(numLonsI - 1)
    dyI = (yValues(numLats) - yValues(1))/(numLatsI - 1)

XValuesI = (/ ( (xValues(1) + dxI*(j-1)), j=1,numLonsI ) /)

yValuesI = (/ ( (yValues(1) + dyI*(j-1)), j=1,numLatsI ) /)

!Index to know which points of the original data use for interpolation
allocate(indi(numLonsI), stat = status)
allocate(indj(numLatsI), stat = status)

do k = 1, numLonsI
indi(k) = dxI/dx*(k-1)!rounded to integer. the order is kinda important
enddo
do j = 1, numLatsI
indj(j) = dyI/dy*(j-1)!rounded to integer. the order is kinda important
enddo

!Keep the left, right, bottom and top boundaries as they were
ZValuesI(1, 1) = ZValues(1, 1)
ZValuesI(1, numLatsI) = ZValues(1, numLats)
ZValuesI(numLonsI, 1) = ZValues(numLons, 1)
ZValuesI(numLonsI, numLatsI) = ZValues(numLons, numLats)

!     Careful! no conversion to reals. Rounding errors are all over the place
     do j = 2, numLatsI-1
          weight1 = (yValuesI(j) - yValues(indj(j)))/dy
          weight2 = (yValues(indj(j)+1) - yValuesI(j))/dy

          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 1', weight1, weight2);
                stop
          end if

          ZValuesI(1, j) = ZValues(1, indj(j) + 1)*weight1 + &
                               ZValues(1, indj(j))*weight2
     enddo
     do j = 2, numLatsI-1
          weight1 = (yValuesI(j) - yValues(indj(j)))/dy
          weight2 = (yValues(indj(j)+1) - yValuesI(j))/dy
          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 2', weight1, weight2);
                stop
          end if

          ZValuesI(numLonsI, j) = ZValues(numLons, indj(j) + 1)*weight1 + &
                               ZValues(numLons, indj(j))*weight2
     enddo
     do k = 2, numLonsI-1
          weight1 = (xValuesI(k) - xValues(indi(k)))/dx
          weight2 = (xValues(indi(k)+1) - xValuesI(k))/dx
          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 3', weight1, weight2);
                stop
          end if

          ZValuesI(k, 1) = ZValues(indi(k) + 1, 1)*weight1 + &
                               ZValues(indi(k), 1)*weight2
     enddo
     do k = 2, numLonsI-1
          weight1 = (xValuesI(k) - xValues(indi(k)))/dx
          weight2 = (xValues(indi(k)+1) - xValuesI(k))/dx
          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 4', weight1, weight2);
                stop
          end if

          ZValuesI(k, numLatsI) = ZValues(indi(k) + 1, numLats)*weight1 + &
                               ZValues(indi(k), numLats)*weight2
     enddo

    do k = 2, numLonsI-1
         do j = 2, numLatsI-1
         !     First horizontal averaging
          weight1 = (xValuesI(k) - xValues(indi(k)))/dx
          weight2 = (xValues(indi(k)+1) - xValuesI(k))/dx
          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging 5', weight1, weight2);
                stop
          end if

          ztmp1 = ZValues(indi(k) + 1, indj(j))*weight1 + &
                              ZValues(indi(k), indj(j))*weight2
          ztmp2 = ZValues(indi(k) + 1, indj(j)+1)*weight1 + &
                              ZValues(indi(k), indj(j)+1)*weight2

         !     Then vertical
          weight1 = (yValuesI(j) - yValues(indj(j)))/dy
          weight2 = (yValues(indj(j)+1) - yValuesI(j))/dy
          if (weight1<0 - eps_tmp .or. weight1>1 + eps_tmp .or. weight2<0 - eps_tmp .or. weight2>1 + eps_tmp )  then
                call handle_av_err('Incorrect averaging', weight1, weight2);
                stop
          end if

          ZValuesI(k, j) = ztmp2*weight1 + ztmp1*weight2
         enddo
    enddo


    end subroutine
!********************************************************************
!********************************************************************
!     Error Management
!********************************************************************
!********************************************************************

!     subroutine handle_nc_err(err_code, section)
!
!     integer, intent(in)                    ::     err_code
!     character (len = *), intent(in) :: section
!
!     if (err_code /= nf90_NoErr) then
!           write(*,*) 'Something went wrong in section ', section
!    end if
!
!     end subroutine
!
!
!     subroutine handle_av_err(message, weight1, weight2)
!
!     character (len = *), intent(in) :: message
!     real(dp)     ::     weight1, weight2
!
!           print *,  message, ", using weights: ", weight1, weight2
!
!     end subroutine
!
!********************************************************************

!********************************************************************
!-------------------- OLD VERSION, OBSOLETE -------------------------
!********************************************************************


end

