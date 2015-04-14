module SHRotate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This module contains an interface block defining all the routines
!	used in the archive SHTOOLS. These are necessary in order to use
!	implicitly shaped arrays with most subroutines.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use my_trigs


	integer, parameter ::	CSPHASE_DEFAULT = 1	! The default is to EXCLUDE the 
							! CONDON-SHORTLEY phase of (-1)^m
							! in front of the Legendre functions.
							! To use this phase function, set
							! CSPHASE_DEFAULT = -1

contains

subroutine SHRotateRealCoef(cilmrot, cilm, lmax, x, dj)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will rotate a set of real spherical harmonic coefficients 
!	corresponding to the angles listed in the input array x.
!
!	The rotation of a coordinate system or body can be viewed in two complementary ways
!	involving three successive rotations. Both methods have the same initial and final 
!	configurations, and the angles listed in both schemes are the same.
!
!	Scheme A: 	(I) Rotation about the z axis by alpha.
!			(II) Rotation about the new y axis by beta.
!			(III) Rotation about the new z axis by gamma.
!
!	Scheme B:	(I) Rotation about the z axis by gamma.
!			(II) Rotation about the initial y axis by beta.
!			(III) Rotation about the initial z axis by alpha.
!
!	The rotations can further be viewed either as either a rotation of the coordinate system
!	or the physical body.
!
!	1. Rotation of the coordinate system without rotation of the physical body,
!		use x(alpha, beta, gamma).
!
!	2. Rotation of the physical body without rotation of the coordinate system,
!		use x(-gamma, -beta, -alpha).
!
!	To perform the inverse trasform of x(alpha, beta, gamma), use x(-gamma, -beta, -alpha).
!
!	This routine uses the "y-convention" were rotations are about the y-axis instead 
!	of the x-axis.
!
!	Calling Parameters
!		IN
!			cilm		Real "geodesy" normalized spherical harmonic coefficients
!					with dimension (2, lmax+1, lmax+1).
!			x		Array or rotation angles.
!			lmax		Maximum spherical harmonic degree.
!		OUT
!			cilmrot		Rotated real "geodesy" normalized spherical harmonic 
!					coefficients with dimension (2, lmax+1, lmax+1).
!		OPTIONAL
!			dj		Rotation matrix with dimension (lmax+1, lmax+1, lmax+1).
!
!	Note: Before using this routine, I would verify that the input euler angles and signs 
!		give the expected results. Some people define the angle beta as a rotation with 
!		respect to the x axis.
!
!	Written by Mark Wieczorek (November 2003)
!	June 4 2006. The routine now checks to see what the default Condon-Shortley phase convention is
!	before converting to complex form for use in the rotation routines.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	implicit none
	
	real*8, intent(in) ::	cilm(:,:,:), x(:), dj(:,:,:)
	real*8, intent(out) ::	cilmrot(:,:,:)
	integer, intent(in) ::	lmax
	integer ::	astat(3)
	real*8, allocatable ::	ccilm(:,:,:), cof(:,:), rcof(:,:)
	
	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- SHRotateRealCoef"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is", lmax
		print*, "Input array is dimensioned ", size(cilm(:,1,1)),  size(cilm(1,:,1)), size(cilm(1,1,:)) 
		stop
	elseif (size(cilmrot(:,1,1)) < 2 .or. size(cilmrot(1,:,1)) < lmax+1 .or. size(cilmrot(1,1,:)) < lmax+1) then
		print*, "Error --- SHRotateRealCoef"
		print*, "CILMROT must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is", lmax
		print*, "Input array is dimensioned ", size(cilmrot(:,1,1)),  size(cilmrot(1,:,1)), size(cilmrot(1,1,:)) 
		stop
	elseif (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax+1 .or. size(dj(1,1,:)) < lmax+1) then
		print*, "Error --- SHRotateRealCoef"
		print*, "DJ must be dimensioned as (LMAX+1, LMAX+1, LMAX+1) where LMAX is", lmax
		print*, "Input array is dimensioned ", size(dj(:,1,1)),  size(dj(1,:,1)), size(dj(1,1,:)) 
		stop
	elseif (size(x) < 3) then
		print*, "Error --- SHRotateRealCoef"
		print*, "X must be dimensioned as (3)"
		print*, "Input array is dimensioned ", size(x)
		stop
	endif
	
	allocate(ccilm(2,lmax+1,lmax+1), stat = astat(1))
	allocate(cof(2,(lmax+1)*(lmax+2)/2), stat = astat(2))
	allocate(rcof(2,(lmax+1)*(lmax+2)/2), stat = astat(3))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
		print*, "Error --- SHRotateRealCoef"
		print*, "Problem allocating arrays CCILM, COF, and RCOF", &
			astat(1), astat(2), astat(3)
		stop
	endif

	ccilm = 0.0d0
	cof = 0.0d0
	rcof = 0.0d0
	cilmrot = 0.0d0
	
	if (CSPHASE_DEFAULT == 1) then
		call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=1)		! Convert geodesy coefficients to Varshalovich et al. complex form
	else
		call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=0)	
	endif
	
	call SHcilmtocindex(ccilm, cof, lmax)			! Re-order complex coefficients to form a vector
							
	call SHRotateCoef(x, cof, rcof, dj, lmax)		! Rotate complex re-ordered coefficients
			
	call SHcindextocilm(rcof, ccilm, lmax)			! Convert ordered coefficients back to an array
	
	if (CSPHASE_DEFAULT == 1) then
		call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=1)	! Convert Varshalovich et al complex coefficients back to geodesy form
	else
		call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=0)
	endif
	
	deallocate(ccilm)
	deallocate(cof)
	deallocate(rcof)
	
end subroutine SHRotateRealCoef

subroutine SHRotateCoef(x, cof, rcof, dj, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will rotate a set of spherical harmonic coefficients (see 
!	documentation in SHRotateRealCoef for a description of how rotation angles 
!	are defined.) All angles are measured in radians, and spherical harmonic
!	coefficients are fully normalized as in Edmonds, 1960. Note that euler 
!	angles are different in Goldstein and Edmonds.
!
!	This routine uses the "y-convention" where rotations are about the y axis
!	instead of the x axis.
!
!	djpi2 must be called before using this routine.
!
!	Calling Parameters
!		IN
!			x	Array of dimension 3 containing the three Euler angles.
!			dj	Roation matrix with dimension (lmax+1, lmax+1, lmax+1).
!			cof	Indexed spherical harmonic coefficients with dimensions
!				(2, (lmax+1)*(lmax+2)/2).
!			lmax	Maximum spherical harmonic degree.
!		OUT
!			rcof	Indexed rotated spherical harmonic coefficients with 
!				dimensions (2, (lmax+1)*(lmax+2)/2).
!
!	Dependencies:	None
!
!	History:
!		1.	Based on routine from Guy Masters (July16, 1993)
!		2.	Modified by Mark Simons (July 25, 1993)
!		3.	Turned into readable f95 code by Mark Wieczorek (August, 2003)
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: 	cof(:,:), dj(:,:,:), x(3)
	real*8, intent(out) ::	rcof(:,:)
	integer, intent(in) :: lmax
	real*8 ::	sum(2), temp(2,lmax+1), temp2(2,lmax+1), cgam(lmax+1), sgam(lmax+1), &
			calf(lmax+1), salf(lmax+1), cbet(lmax+1), sbet(lmax+1), pi2, alpha, &
			beta, gamma
	integer ::	ind, lp1, l, mp1, jp1, isgn, ii, indx
	
	if (size(cof(:,1)) < 2 .or. size(cof(1,:)) < (lmax+1)*(lmax+2)/2) then
		print*, "Error --- SHRotateCoef"
		print*, "COEF must be dimensioned (2, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(cof(:,1)),  size(cof(1,:))
		stop
	elseif (size(rcof(:,1)) < 2 .or. size(rcof(1,:)) < (lmax+1)*(lmax+2)/2) then
		print*, "Error --- SHRotateCoef"
		print*, "RCOEF must be dimensioned (2, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(rcof(:,1)),  size(rcof(1,:))
		stop
	elseif (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax+1 .or. size(dj(1,1,:)) < lmax+1 ) then
		print*, "Error --- SHRotateCoef"
		print*, "DJ must be dimensioned (LMAX+1, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(dj(:,1,1)), size(dj(1,:,1)), size(dj(1,1,:))
		stop
	endif

	pi2 = 1.570796326794895d0

	alpha = x(1)
	beta  = x(2)
	gamma = x(3)
	
	alpha = alpha-pi2
	gamma = gamma+pi2
	beta = -beta
	
	ind = 0
	
	! Loop over degrees
	
	do lp1 = 1, lmax+1 
		l = lp1-1
		cbet(lp1) = cos(l*beta)
		sbet(lp1) = sin(l*beta)
		cgam(lp1) = cos(l*gamma)
		sgam(lp1) = sin(l*gamma)
		calf(lp1) = cos(l*alpha)
		salf(lp1) = sin(l*alpha)

		! Alpha rotation
		
		do mp1 = 1, lp1
			indx = ind+mp1
			temp(1,mp1) = cof(1,indx) * calf(mp1) - cof(2,indx) * salf(mp1)
			temp(2,mp1) = cof(2,indx) * calf(mp1) + cof(1,indx) * salf(mp1)
		enddo
		
		! B rotation and beta rotation
		
		do jp1 = 1, lp1
			sum(1) = dj(jp1,1,lp1)*temp(1,1)
			sum(2) = 0.d0
			isgn = 1-2*mod((lp1-jp1),2)
			
			do mp1 = 2, lp1
				isgn = -isgn
				ii = (3-isgn)/2
				sum(ii) = sum(ii)+2.d0*dj(jp1,mp1,lp1)*temp(ii,mp1)
			enddo
			
			temp2(1,jp1) = sum(1)*cbet(jp1)-sum(2)*sbet(jp1)
			temp2(2,jp1) = sum(2)*cbet(jp1)+sum(1)*sbet(jp1)
		
		enddo

		! Inverse B rotation and gamma rotation
		
		do jp1 = 1, lp1
		
			sum(1) = dj(1,jp1,lp1)*temp2(1,1)
			sum(2) = 0.d0
			isgn = 1-2*mod((lp1-jp1),2)
			
			do mp1 = 2,lp1
				isgn = -isgn
				ii = (3-isgn)/2
				sum(ii) = sum(ii)+2.d0*dj(mp1,jp1,lp1)*temp2(ii,mp1)
			enddo
			
			indx = ind+jp1
			rcof(1,indx) = sum(1)*cgam(jp1)-sum(2)*sgam(jp1)
			rcof(2,indx) = sum(2)*cgam(jp1)+sum(1)*sgam(jp1)
			
		enddo
		
		ind = ind+lp1
		
	enddo
      
end subroutine SHRotateCoef

subroutine SHrtoc(rcilm, ccilm, degmax, convention, switchcs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert the real "geodesy 4-pi" spherical harmonic coefficients into 
!	complex form with either a 4 pi or unit normalization. The complex 
!	coefficients are only calculated for positive m's (the negative m's are given by 
!
!		c_l,m* = (-1)^m c_l,-m
!
!	This sign convention is designed to be used with the spherical haromic rotation routines
!	taken from Mark Simons' code. If degmax is not specified, then the maximum degree of the 
!	conversion is taken from the size of the input arrays.
!
!
!	Calling Parameters
!		IN
!			rcilm		Real "geodesy" spherical harmonic coefficients with dimensions
!					(2, lmax+1, lmax+1).
!		OUT
!			ccilm		Complex unity-normalized spherical harmonic coefficients, dimensioned
!					as (2, lmax+1, lmax+1). The first index corresponds to the real and
!					complex coefficients, respectively.
!		OPTIONAL
!			degmax		Maximum degree of conversion to be performed.
!			convention	1=output 4-pi normalized coefficients
!					2=output  Varshalovich et al. normalized coefficients
!			switchcs	If 1, Change between different Condon-Shortley phase convenctions.
!					If 0, use consistent phase convention.
!
!		Dependencies:	None
!			
!
!	Written by Mark Wieczorek (August 2003)
!	June 4 2006: Added optional argument to change between Condon-Shortley phase conventions.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: 		rcilm(:,:,:)
	real*8, intent(out) ::	 	ccilm(:,:,:) 
	integer, intent(in), optional ::	degmax, convention, switchcs
	integer :: 			lmax, l, m
	real*8 ::			pi
	
	
	if (present(switchcs)) then
		if (switchcs /= 1 .and. switchcs /=0) then
			print*, "Error --- SHrtoc"
			print*, "switchcs must be equal to either 0 (keep same convention) of 1 (change Condon-Shortley phase)"
			print*, "Input value is ", switchcs
			stop
		endif
	endif
	
	if (present(convention) ) then
		if (convention /=1 .and. convention /=2) then
			print*, "Error --- SHrtoc"
			print*, "CONVENTION must be 1 or 2."
			print*, "Input valuse is ", convention
			stop
		endif
	endif

	if (present(degmax)) then
		lmax = degmax
		if (size(rcilm(:,1,1)) < 2 .or. size(rcilm(1,:,1)) < lmax +1 .or. size(rcilm(1,1,:)) < lmax +1) then
			print*, "Error --- SHrtoc"
			print*, "RCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned as ", size(rcilm(:,1,1)), size(rcilm(1,:,1)),  size(rcilm(1,1,:))
			stop
		elseif (size(ccilm(:,1,1)) < 2 .or. size(ccilm(1,:,1)) < lmax +1 .or. size(ccilm(1,1,:)) < lmax +1) then
			print*, "Error --- SHrtoc"
			print*, "CCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned as ", size(ccilm(:,1,1)), size(ccilm(1,:,1)),  size(ccilm(1,1,:))
			stop

		endif
	else
		if (size(rcilm(:,1,1)) < 2) then
			print*, "Error --- SHrtoc"
			print*, "RCILM must be dimensioned as (2,*,*)."
			print*, "Input array is dimensioned as ",  size(rcilm(:,1,1)), size(rcilm(1,:,1)), size(rcilm(1,1,:))
			stop
		elseif (size(ccilm(:,1,1)) < 2) then
			print*, "Error --- SHrtoc"
			print*, "CCILM must be dimensioned as (2,*,*)."
			print*, "Input array is dimensioned as ",  size(ccilm(:,1,1)), size(ccilm(1,:,1)), size(ccilm(1,1,:))
			stop
		endif
		
		lmax = min(size(rcilm(1,1,:)) -1, size(ccilm(1,1,:)) -1, size(rcilm(1,:,1)) -1, size(ccilm(1,:,1)) -1)
					
	endif

	pi = acos(-1.0d0)
	ccilm= 0.0d0

		
	do l = 0, lmax, 1
	
		if (present(convention) ) then
			if (convention == 2) then
				ccilm(1,l+1, 1) = sqrt(4.0d0*pi)*rcilm(1,l+1,1)
				ccilm(2,l+1, 1) = 0.0d0
			
				do m=1, l, 1
					if (present(switchcs)) then
						if (switchcs == 1) then
							ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) * (-1.0d0)**m
							ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1) * (-1.0d0)**m
						else
							ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) 
							ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1)
						endif
					else
						ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) 
						ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1)
					endif
				enddo
			endif
			
		else
			ccilm(1,l+1, 1) = rcilm(1,l+1,1)
			ccilm(2,l+1, 1) = 0.0d0
			
			do m=1, l, 1
				if (present(switchcs)) then
					if (switchcs == 1) then
						ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1)/sqrt(2.0d0) * (-1.0d0)**m
						ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1)/sqrt(2.0d0) * (-1.0d0)**m
					else
						ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1)/sqrt(2.0d0)
						ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1)/sqrt(2.0d0)
					endif
				else
					ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1)/sqrt(2.0d0)
					ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1)/sqrt(2.0d0)
				endif
			enddo
		endif
				
	enddo
	
end subroutine SHrtoc


subroutine SHctor(ccilm, rcilm, degmax, convention, switchcs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert either "geodesy 4-pi" or "Varshalovich et al." complex spherical 
!	harmonic coefficients into real "geodesy 4-pi" spherical harmonic coefficients.
!
!	If degmax is not specified, then the maximum degree of the 
!	conversion is taken from the size of the input arrays.
!
!
!	Calling Parameters
!		IN
!			ccilm		Complex unity-normalized spherical harmonic coefficients, dimensioned
!					as (2, lmax+1, lmax+1). The first index corresponds to the real and
!					complex coefficients, respectively.
!		OUT
!			rcilm		Real "geodesy" spherical harmonic coefficients with dimensions
!					(2, lmax+1, lmax+1).
!		OPTIONAL
!			degmax		Maximum degree of conversion to be performed.
!			convention	1=input coefficients are 4-pi normalized
!					2=input coefficients are Varshalovich et al. normalized.
!			switchcs	If 1, Change between different Condon-Shortley phase convenctions.
!					If 0, use consistent phase convention.
!
!		Dependencies:	None
!
!	Written by Mark Wieczorek 2004.
!	June 4 2006: Added option to change between Condon-Shortley phase conventions.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: ccilm(:,:,:)
	real*8, intent(out) :: rcilm(:,:,:)
	integer, intent(in), optional ::	degmax, convention, switchcs
	integer ::	lmax, l, m
	real*8 ::	pi
	
	
	if (present(switchcs)) then
		if (switchcs /= 1 .and. switchcs /=0) then
			print*, "Error --- SHrtoc"
			print*, "switchcs must be equal to either 0 (keep same convention) of 1 (change Condon-Shortley phase)"
			print*, "Input value is ", switchcs
			stop
		endif
	endif
	
	if (present(convention) ) then
		if (convention /=1 .and. convention /=2) then
			print*, "Error --- SHrtoc"
			print*, "CONVENTION must be 1 or 2."
			print*, "Input valuse is ", convention
			stop
		endif
	endif

	if (present(degmax)) then
		lmax = degmax
		if (size(rcilm(:,1,1)) < 2 .or. size(rcilm(1,:,1)) < lmax +1 .or. size(rcilm(1,1,:)) < lmax +1) then
			print*, "Error --- SHrtoc"
			print*, "RCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned as ", size(rcilm(:,1,1)), size(rcilm(1,:,1)),  size(rcilm(1,1,:))
			stop
		elseif (size(ccilm(:,1,1)) < 2 .or. size(ccilm(1,:,1)) < lmax +1 .or. size(ccilm(1,1,:)) < lmax +1) then
			print*, "Error --- SHrtoc"
			print*, "CCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned as ", size(ccilm(:,1,1)), size(ccilm(1,:,1)),  size(ccilm(1,1,:))
			stop

		endif
	else
		if (size(rcilm(:,1,1)) < 2) then
			print*, "Error --- SHrtoc"
			print*, "RCILM must be dimensioned as (2,*,*)."
			print*, "Input array is dimensioned as ",  size(rcilm(:,1,1)), size(rcilm(1,:,1)), size(rcilm(1,1,:))
			stop
		elseif (size(ccilm(:,1,1)) < 2) then
			print*, "Error --- SHrtoc"
			print*, "CCILM must be dimensioned as (2,*,*)."
			print*, "Input array is dimensioned as ",  size(ccilm(:,1,1)), size(ccilm(1,:,1)), size(ccilm(1,1,:))
			stop
		endif
		
		lmax = min(size(rcilm(1,1,:)) -1, size(ccilm(1,1,:)) -1, size(rcilm(1,:,1)) -1, size(ccilm(1,:,1)) -1)
					
	endif
	
	
	pi = acos(-1.0d0)
	rcilm = 0.0d0

		
	do l = 0, lmax, 1
	
		if (present(convention) ) then
			if (convention == 2) then
				rcilm(1,l+1,1) = ccilm(1,l+1,1)/sqrt(4.0d0*pi)
				rcilm(2, l+1, 1) = 0.0d0
			
				do m=1, l, 1
					if (present(switchcs)) then
						if (switchcs == 1) then
							rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi) * (-1.0d0)**m
							rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) /sqrt(2.0d0*pi) * (-1.0d0)**m
						else
							rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi)
							rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) /sqrt(2.0d0*pi)
						endif
					else
						rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi)
						rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) /sqrt(2.0d0*pi) 
					endif
				enddo
			endif
			
		else
		
			rcilm(1,l+1,1) = ccilm(1,l+1,1)
			rcilm(2, l+1, 1) = 0.0d0
			
			do m=1, l, 1
				if (present(switchcs)) then
					if(switchcs == 1) then
						rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1) * (-1.0d0)**m
						rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1) * (-1.0d0)**m
					else
						rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1)
						rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1)
					endif
				else
					rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1)
					rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1)
				endif
			enddo
		
		endif
	enddo
	
end subroutine SHctor


subroutine SHCilmToCindex(cilm, cindex, degmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will convert a 3D matrix of spherical harmonics indexed as (i, l+1, m+1)
!	into a 2D matrix that is indexed as (i, index) where index = l(l+1)/2+m+1.
!
!	Calling Parameters:
!		IN
!			cilm	Array of spherical harmonic coefficients with dimensions 
!				(2, lmax+1, lmax+1).
!		OUT
!			cindex	Array of indexed spherical harmonic coefficnets with dimensions
!				(2, (lmax+1)*(lmax+2)/2).
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: cilm(:,:,:)
	real*8, intent(out) :: cindex(:,:)
	integer, intent(in), optional ::	degmax
	integer :: lmax, l, m, index
	
	if (present(degmax)) then
		lmax = degmax
		if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax + 1 .or. size(cilm(1,1,:)) < lmax+1 ) then
			print*, "Error --- SHcilmtocindex"
			print*, "CILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
			stop
		elseif (size(cindex(:,1)) < 2 .or. size(cindex(1,:)) < (lmax+1)*(lmax+2)/2) then
			print*, "Error --- SHcilmtocindex"
			print*, "CINDEX must be dimensioned as (2, (DEGMAX+1)*(DEGMAX+2)/2) where DEGMAX is ", degmax
			print*, "Input array is dimensioned ", size(cindex(:,1)), size(cindex(1,:))
			stop
		endif
	else
		lmax = min(size(cilm(1,1,:)) - 1, size(cilm(1,:,1)) - 1)
		if (size(cilm(:,1,1)) < 2) then
			print*, "Error --- SHcilmtocindex"
			print*, "CILM must be dimensioned as (2, *, *)."
			print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
			stop
		elseif (size(cindex(:,1)) < 2 .or. size(cindex(1,:)) <  (lmax+1)*(lmax+2)/2) then
			print*, "Error --- SHcilmtocindex"
			print*, "CINDEX must be dimensioned as (2, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
			stop
		endif
	endif
		
	do l=0, lmax
		do m=0, l
			index = l*(l+1)/2+m+1
			cindex(1, index) = cilm(1,l+1,m+1)
			cindex(2, index) = cilm(2,l+1,m+1)
		enddo
	enddo

end subroutine SHCilmToCindex
		
		
subroutine SHCindexToCilm(cindex, cilm, degmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will convert a 2D matrix of spherical harmonics indexed as (i, index)
!	into a 2D matrix that is index as (i, l+1, m+1) where index = l(l+1)/2+m+1.
!
!	Calling Parameters:
!		IN
!			cindex	Array of indexed spherical harmonic coefficnets with dimensions
!				(2, (lmax+1)*(lmax+2)/2).
!		OUT
!			cilm	Array of spherical harmonic coefficients with dimensions 
!				(2, lmax+1, lmax+1).
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(out) :: cilm(:,:,:)
	real*8, intent(in) ::	cindex(:,:)
	integer, intent(in), optional ::	degmax
	integer :: lmax, l, m, index, n
	
	n = size(cindex(1,:))

	if (present(degmax)) then
		lmax = degmax
		if (lmax > nint((-3.0d0 + sqrt(9.0d0 + 8.0d0*dble(n-2)) )/2.0d0) ) then
			print*, "Error - SHcindextocilm"
			print*, "The output spherical harmonic degree DEGMAX is larger than the input coefficients."
			print*, "Input value of DEGMAX ", degmax
			print*, "Maximum spherical harmonic degree of CINDEX ", nint((-3.0d0 + sqrt(9.0d0 + 8.0d0*dble(n-2)) )/2.0d0)
			stop
		elseif (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax + 1 .or. size(cilm(1,1,:)) < lmax+1 ) then
			print*, "Error --- SHcindextocilm"
			print*, "CILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
			stop
		endif
	else	
		lmax = nint((-3.0d0 + sqrt(9.0d0 + 8.0d0*dble(n-2)) )/2.0d0)
		
		if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax + 1 .or. size(cilm(1,1,:)) < lmax+1 ) then
			print*, "Error --- SHcindextocilm"
			print*, "CILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) where DEGMAX is ", degmax
			print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
			stop
		endif
	endif
		
	do l=0, lmax
		do m=0, l
			index = l*(l+1)/2+m+1
			cilm(1,l+1,m+1) = cindex(1, index)
			cilm(2,l+1,m+1) = cindex(2, index)	
		enddo
	enddo	

end subroutine SHCindexToCilm

subroutine djpi2(dj, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine computes 
!	  j
!	d     (pi/2)
!	  m N
!
!	for all posible values of m (>=0) and N for 0<l<lmax.
!	The output array corresponds to dj(N,m,j). (I'm not positive
!	about the N and m ordering).
!
!	Note that this algorithm is easy to modify to work on a
!	single l (just need a loop to compute f1) since there are
!	no l recursions involved.
!
!	Calling Parameters
!		IN
!			lmax	Maximum spherical harmonic degree to be computed.
!		OUT
!			dj	Rotation matrix with dimensions (lmax+1, lmax+1, lmax+1).
!
!	Dependencies:	None
!
!	History
!
!		1. Based on routine from Guy Masters (July 16, 1993)
!		2. Modified by Mark Simons (July 19, 1993)
!		3. Turned into readable f95 code by Mark Wieczorek (August, 2003).
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	dj(:,:,:)
	integer ::		i, l, n, lp1, j, m, isn, np
	real*8 ::		f((lmax+1)*8), f1, f2, g1, g2, en2, fl2p1
	
	if (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax + 1 .or. size(dj(1,1,:)) < lmax +1) then
		print*, "Error --- djpi2"
		print*, "DJ must be dimensioned (LMAX+1, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(dj(:,1,1)), size(dj(1,:,1)), size(dj(1,1,:))
		stop
	endif
	
	dj = 0.0d0
	
	dj(1,1,1) = 1.d0
	dj(1,1,2) = 0.d0
	dj(1,2,2) = -1.0d0/sqrt(2.0d0)
	dj(2,1,2) = -dj(1,2,2)
	dj(2,2,2) = 0.5d0
	f1        = 0.5d0

	do l = 2, lmax

		lp1   = l+1
		fl2p1 = l+lp1
		
		do i = 1, l
			f(i) = dsqrt(i*(fl2p1-i))
        	enddo
        	
        	f1 = f1*(l+l-1.0d0)/(l+l)

		! Do N = 0 terms
		
		dj(lp1,1,lp1) = -dsqrt(f1)
		dj(l,1,lp1)   = 0.0d0
		
		do i = 2, l
			j = lp1-i
			dj(j,1,lp1) = -f(i-1)*dj(j+2,1,lp1)/f(i)
		enddo

		! Do positive N terms (bottom triangle)
		
		f2 = f1
		g1 = l
		g2 = lp1
		
		do n = 1, l
			np = n+1
			en2 = n+n
			g1 = g1 + 1.0d0
			g2 = g2 - 1.0d0
			f2 = f2*g2/g1
			dj(lp1,np,lp1) = -dsqrt(f2)
			dj(l,np,lp1) = dj(lp1,np,lp1)*en2/f(1)

			do i = 2, l-n
				j = lp1-i
				dj(j,np,lp1) = ( en2*dj(j+1,np,lp1) &
					- f(i-1)*dj(j+2,np,lp1) ) / f(i)
			enddo
			
		enddo

		! Fill upper triangle and fix signs

		do j = 1, l
			do m = j,l
				dj(j,m+1,lp1) = dj(m+1,j,lp1)
			enddo
		enddo

		isn = 1 + mod(l,2)
		
		do np = 1, lp1
			do i = isn, lp1, 2
				dj(i,np,lp1) = -dj(i,np,lp1)
			enddo
		enddo
		
	enddo
	
end subroutine djpi2


subroutine SHRead(filename, cilm, lmax, skip, header, error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will open an ascii file, and read in all of the 
!	spherical harmonic coefficients. If the option "header"
!	is specified, then the first Length(header) records of the 
!	spherical harmonic file will be retuned in the array header.
!	
!	The file will be read until the end of the file is encountered
!	or until the maximum length of cilm (as dimensioned in the calling 
!	program) is reached.
!
!	Calling Parameters:
!		IN
!			filename	Character name of the ascii file.
!		OUT
!			cilm		Spherical harmonic coeficients with dimensions
!					(2, lmax+1, lmax+1).
!			lmax		Maximum spherical harmonic degree of cilm.
!		OPTIONAL
!			header		Array of the first length(header) data records
!					in the file. Called as header = header, or 
!					header = header(1:number_of_header_records).
!			skip		Number of lines to skip
!			error 		Error of the spherical harmonic coefficients, assumed 
!					to in the format (l, m, c1lm, c2lm, error1lm, error2lm).
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek (September 2003)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	character(*), intent(in)::		filename
	integer, intent(out) ::			lmax
	real*8, intent(out) ::			cilm(:,:,:)
	real*8, intent(out), optional ::	header(:), error(:,:,:)
	integer, intent(in), optional ::	skip
	integer ::				l, m, stat, ll, mm, lmax2, lstart, headlen, fu
	
	lmax = 0
	cilm = 0.0d0
	fu = 101
	
	if (size(cilm(:,1,1)) < 2 ) then
		print*, "Error --- SHRead"
		print*, "CILM must be dimensioned (2, *, *)."
		print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	lmax2 = min(size(cilm(1,1,:) ) - 1, size(cilm(1,:,1) ) - 1)
	
	open(fu, file=filename, status="old")
		
	if (present(header)) headlen = size(header)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	SKip lines and read header information
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (present(skip) ) then
		do l=1, skip, 1
			read(fu,*, iostat=stat)
			if (stat /= 0 ) then
				print*, "Error --- SHRead"
				print*, "Problem skipping first lines of ", filename
				print*, "Line number = ", l
				print*, "Number of lines to skip = ", skip
				stop
			endif
		enddo
	endif
	
	if (present(header) ) then
		read(fu,*, iostat=stat) (header(l), l=1, headlen)
		if (stat /= 0 ) then
			print*, "Error --- SHRead"
			print*, "Problem reading header line ", filename 
			stop
		endif
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine first l value
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	read(fu,*, iostat=stat) ll
	
	if (stat /= 0 ) then
		print*, "SHRead --- Error "
		print*, "Problem reading first line of ", filename
		stop
	endif

	lstart = ll
	
	rewind(fu)
	
	if ( present(skip) ) then
		do l=1, skip, 1
			read(fu,*, iostat=stat)
		enddo
	endif
	
	if ( present(header) ) then
		read(fu,*, iostat=stat) (header(l), l=1, headlen)
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Read coefficients
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do l=lstart, lmax2, 1
       		do m=0,l
       			if ( present(error) ) then
       				read(fu,*,iostat=stat) ll, mm, cilm(1,l+1,m+1), cilm(2,l+1,m+1), &
       					error(1,l+1,m+1), error(2,l+1,m+1)
       			else
       				if (m==0) then
       					read(fu,*,iostat=stat) ll, mm, cilm(1,l+1,m+1)
       				else
       					read(fu,*,iostat=stat) ll, mm, cilm(1,l+1,m+1), cilm(2,l+1,m+1)
       				endif
       			endif
             		if (stat < 0) then
       				exit
       			elseif (stat > 0) then
       				print*, "SHRead --- Error "
				print*, "Problem reading file ", filename
				stop
			elseif (ll /=l .or. mm /=m) then
				print*, "SHRead --- Error "
				print*, "Problem reading file ", filename
				print*, "Expected indices (l,m) = ", l, m
				print*, "Read indices (l,m) = ", ll, mm
				stop
			endif
		enddo
		
		if (stat < 0) exit
	
		lmax = l	
		
	enddo
	
	close(fu)
		
end subroutine SHRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SHRotateRealCoef_2(cilmrot, cilm, lmax, x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will rotate a set of real spherical harmonic coefficients
!	corresponding to the angles listed in the input array x.
!	SAME AS SHRotateRealCoef, EXCEPT THAT SIZE(LMAX^3) DJ IS NOT STORED.
!
!	Calling Parameters
!		IN
!			cilm		Real "geodesy" normalized spherical harmonic coefficients
!					with dimension (2, lmax+1, lmax+1).
!			x		Array or rotation angles.
!			lmax		Maximum spherical harmonic degree.
!		OUT
!			cilmrot		Rotated real "geodesy" normalized spherical harmonic
!					coefficients with dimension (2, lmax+1, lmax+1).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	implicit none

	real*8, intent(in) ::	cilm(:,:,:), x(:)
	real*8, intent(out) ::	cilmrot(:,:,:)
	integer, intent(in) ::	lmax
	integer ::	astat(3)
	real*8, allocatable ::	ccilm(:,:,:), cof(:,:), rcof(:,:)

	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- SHRotateRealCoef"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is", lmax
		print*, "Input array is dimensioned ", size(cilm(:,1,1)),  size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	elseif (size(cilmrot(:,1,1)) < 2 .or. size(cilmrot(1,:,1)) < lmax+1 .or. size(cilmrot(1,1,:)) < lmax+1) then
		print*, "Error --- SHRotateRealCoef"
		print*, "CILMROT must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is", lmax
		print*, "Input array is dimensioned ", size(cilmrot(:,1,1)),  size(cilmrot(1,:,1)), size(cilmrot(1,1,:))
		stop
	elseif (size(x) < 3) then
		print*, "Error --- SHRotateRealCoef"
		print*, "X must be dimensioned as (3)"
		print*, "Input array is dimensioned ", size(x)
		stop
	endif

	allocate(ccilm(2,lmax+1,lmax+1), stat = astat(1))
	allocate(cof(2,(lmax+1)*(lmax+2)/2), stat = astat(2))
	allocate(rcof(2,(lmax+1)*(lmax+2)/2), stat = astat(3))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
		print*, "Error --- SHRotateRealCoef"
		print*, "Problem allocating arrays CCILM, COF, and RCOF", &
			astat(1), astat(2), astat(3)
		stop
	endif

	ccilm = 0.0d0
	cof = 0.0d0
	rcof = 0.0d0
	cilmrot = 0.0d0

	if (CSPHASE_DEFAULT == 1) then
		call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=1)		! Convert geodesy coefficients to Varshalovich et al. complex form
	else
		call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=0)
	endif

	call SHcilmtocindex(ccilm, cof, lmax)			! Re-order complex coefficients to form a vector

	call SHRotateCoef_2(x, cof, rcof, lmax)		! Rotate complex re-ordered coefficients

	call SHcindextocilm(rcof, ccilm, lmax)			! Convert ordered coefficients back to an array

	if (CSPHASE_DEFAULT == 1) then
		call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=1)	! Convert Varshalovich et al complex coefficients back to geodesy form
	else
		call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=0)
	endif

	deallocate(ccilm)
	deallocate(cof)
	deallocate(rcof)

end subroutine SHRotateRealCoef_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SHRotateCoef_2(x, cof, rcof, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: 	cof(:,:), x(3)
	real*8, intent(out) ::	rcof(:,:)
	integer, intent(in) :: lmax
	real*8 ::	f1, dj(lmax+1,lmax+1)!, dj_all(lmax+1,lmax+1,lmax+1),
	real*8 ::	sum(2), temp(2,lmax+1), temp2(2,lmax+1), cgam(lmax+1), sgam(lmax+1), &
			calf(lmax+1), salf(lmax+1), cbet(lmax+1), sbet(lmax+1), pi2, alpha, &
			beta, gamma
	integer ::	ind, lp1, l, mp1, jp1, isgn, ii, indx

	if (size(cof(:,1)) < 2 .or. size(cof(1,:)) < (lmax+1)*(lmax+2)/2) then
		print*, "Error --- SHRotateCoef"
		print*, "COEF must be dimensioned (2, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(cof(:,1)),  size(cof(1,:))
		stop
	elseif (size(rcof(:,1)) < 2 .or. size(rcof(1,:)) < (lmax+1)*(lmax+2)/2) then
		print*, "Error --- SHRotateCoef"
		print*, "RCOEF must be dimensioned (2, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(rcof(:,1)),  size(rcof(1,:))
		stop
	endif

!	pi2 = 1.570796326794895d0
	pi2 = pi/2

	alpha = x(1)
	beta  = x(2)
	gamma = x(3)

	alpha = alpha-pi2
	gamma = gamma+pi2
	beta = -beta

	ind = 0
!	call djpi2(dj_all, lmax)
	! Loop over degrees
!	write(*, '(a, i5, a)', advance='no') 'SHRotateCoef_2 loop counter', lmax, ':'
	do lp1 = 1, lmax+1
!		write(*, '(i5)', advance='no') lp1
		l = lp1-1
		call djpi2_2(dj, l, lmax, f1)
!		dj = dj_all(:,:,lp1)

		cbet(lp1) = cos(l*beta)
		sbet(lp1) = sin(l*beta)
		cgam(lp1) = cos(l*gamma)
		sgam(lp1) = sin(l*gamma)
		calf(lp1) = cos(l*alpha)
		salf(lp1) = sin(l*alpha)

		! Alpha rotation

		do mp1 = 1, lp1
			indx = ind+mp1
			temp(1,mp1) = cof(1,indx) * calf(mp1) - cof(2,indx) * salf(mp1)
			temp(2,mp1) = cof(2,indx) * calf(mp1) + cof(1,indx) * salf(mp1)
		enddo

		! B rotation and beta rotation

		do jp1 = 1, lp1
			sum(1) = dj(jp1,1)*temp(1,1)
			sum(2) = 0.d0
			isgn = 1-2*mod((lp1-jp1),2)

			do mp1 = 2, lp1
				isgn = -isgn
				ii = (3-isgn)/2
				sum(ii) = sum(ii)+2.d0*dj(jp1,mp1)*temp(ii,mp1)
			enddo

			temp2(1,jp1) = sum(1)*cbet(jp1)-sum(2)*sbet(jp1)
			temp2(2,jp1) = sum(2)*cbet(jp1)+sum(1)*sbet(jp1)

		enddo

		! Inverse B rotation and gamma rotation

		do jp1 = 1, lp1

			sum(1) = dj(1,jp1)*temp2(1,1)
			sum(2) = 0.d0
			isgn = 1-2*mod((lp1-jp1),2)

			do mp1 = 2,lp1
				isgn = -isgn
				ii = (3-isgn)/2
				sum(ii) = sum(ii)+2.d0*dj(mp1,jp1)*temp2(ii,mp1)
			enddo

			indx = ind+jp1
			rcof(1,indx) = sum(1)*cgam(jp1)-sum(2)*sgam(jp1)
			rcof(2,indx) = sum(2)*cgam(jp1)+sum(1)*sgam(jp1)

		enddo

		ind = ind+lp1

	enddo

!	print *, ""

end subroutine SHRotateCoef_2

subroutine djpi2_2(dj, l, lmax, f1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SAME AS djpi2 BUT DOES NOT REQUIRE O(LMAX^3) RAM.
!!!!!!!!!!!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HAS TO BE CALLED IN THE LOOP IN L:
! 	do l = 0, lmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	implicit none
	integer, intent(in) ::	l, lmax
	real*8, intent(out) ::	dj(:,:)
	integer ::		i, n, lp1, j, m, isn, np
	real*8, intent(inout) :: f1
	real*8 ::		f((lmax+1)*8), f2, g1, g2, en2, fl2p1

	if (size(dj(:,1)) < lmax+1 .or. size(dj(1,:)) < lmax + 1) then
		print*, "Error --- djpi2"
		print*, "DJ must be dimensioned (LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(dj(:,1)), size(dj(1,:))
		stop
	endif

!	INITIALIZE DJ
	dj = 0.0d0

	if (l == 0) then
	dj(1,1) = 1.d0
	f1 = 0.5d0 ! f1 is initialized here and will be changed later
	elseif (l == 1) then
	dj(1,1) = 0.d0
	dj(1,2) = -1.0d0/sqrt(2.0d0)
	dj(2,1) = -dj(1,2)
	dj(2,2) = 0.5d0
	else

!!	do l = 2, lmax

		lp1   = l+1
		fl2p1 = l+lp1

		do i = 1, l
			f(i) = dsqrt(i*(fl2p1-i))
        enddo

        	f1 = f1*(l+l-1.0d0)/(l+l)

		! Do N = 0 terms

		dj(lp1,1) = -dsqrt(f1)
		dj(l,1)   = 0.0d0

		do i = 2, l
			j = lp1-i
			dj(j,1) = -f(i-1)*dj(j+2,1)/f(i)
		enddo

		! Do positive N terms (bottom triangle)

		f2 = f1
		g1 = l
		g2 = lp1

		do n = 1, l
			np = n+1
			en2 = n+n
			g1 = g1 + 1.0d0
			g2 = g2 - 1.0d0
			f2 = f2*g2/g1
			dj(lp1,np) = -dsqrt(f2)
			dj(l,np) = dj(lp1,np)*en2/f(1)

			do i = 2, l-n
				j = lp1-i
				dj(j,np) = ( en2*dj(j+1,np) &
					- f(i-1)*dj(j+2,np) ) / f(i)
			enddo

		enddo

		! Fill upper triangle and fix signs

		do j = 1, l
			do m = j,l
				dj(j,m+1) = dj(m+1,j)
			enddo
		enddo

		isn = 1 + mod(l,2)

		do np = 1, lp1
			do i = isn, lp1, 2
				dj(i,np) = -dj(i,np)
			enddo
		enddo
!!	enddo
	endif

end subroutine djpi2_2
	
end module SHRotate
