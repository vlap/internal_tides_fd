      program test_shallow

      use sphere
c
c  Test program for sphere module - non-linear steady-state geostropic
c  flow in a shallow water model.
c
c  errors should be O(10E-6) or less in single-precision, 10E-7 or less
c  in double precision.
c
c
c     the nonlinear shallow-water equations on the sphere are
c     solved using a spectral method based on the spherical harmonics.
c     the method is described in the paper:
c
c [1] p. n. swarztrauber, spectral transform methods for solving
c     the shallow-water equations on the sphere, p.n. swarztrauber,
c     monthly weather review, vol. 124, no. 4, april 1996, pp. 730-744.
c
c     this program implements test case 3 (steady nonlinear rotated flow)
c     in the paper:
c
c [2] d.l. williamson, j.b. drake, j.j. hack, r. jakob, and
c     p.n. swarztrauber, j. comp. phys., a standard test set
c     for numerical approximations to the shallow-water
c     equations in spherical geometry, j. comp. phys.,
c     vol. 102, no. 1, sept. 1992, pp. 211-224.
c
c definitions:
c
c
c     nlat          number of latitudes 
c     nlon          number of distinct longitudes
c     ntrunc        max wave number
c     omega         rotation rate of earth in radians per second
c     aa            radius of earth in meters
c     pzero         mean height of geopotential
c     uzero         maximum velocity
c     alpha         tilt angle of the rotated grid
c     ncycle        cycle number
c     time          model time in seconds
c     dt            time step
c     lambda        longitude
c     theta         colatitude
c
c   the second dimension of the following two dimensional arrays
c   corresponds to the latitude index with values j=1,...,nlat
c   going from north to south.
c   the second dimension is longitude with values i=1,...,nlon
c   where i=1 corresponds to zero longitude and j=nlon corresponds
c   to 2pi minus 2pi/nlon.
c
c     u(i,j)       east longitudinal velocity component at t=time
c     v(i,j)       latitudinal velocity component at t=time
c     p(i,j)       +pzero = geopotential at t=time
c
c     divg(i,j)    divergence (d/dtheta (cos(theta) v)
c                                          + du/dlambda)/cos(theta)
c     vrtg(i,j)    vorticity  (d/dtheta (cos(theta) u)
c                                          - dv/dlambda)/cos(theta)
c
c     uxact(i,j)   the "exact" longitudinal velocity component
c     vxact(i,j)   the "exact" latitudinal  velocity component
c     uxact(i,j)   the "exact" geopotential
c
      parameter (nlon=128, nlat=(nlon/2)+1, ntrunc=42)
      parameter (nmdim = (ntrunc+1)*(ntrunc+2)/2)
      complex, dimension(nmdim) :: vrtnm,divnm,pnm,scrnm
      complex, dimension(nmdim,3) :: dvrtdtnm,ddivdtnm,dpdtnm
      complex scrm(ntrunc+1,nlat,4)
      real, dimension(nlon,nlat) :: uxact,vxact,pxact,u,v,p,f,
     * ug,vg,pg,vrtg,divg,scrg1,scrg2
      real lambda,lhat,phlt(361),lap(nmdim),gaulats(nlat),weights(nlat)
      integer indxn(nmdim)
      character*3 gridtype

      if (mod(nlat,2) .eq. 0) then
        gridtype = 'GAU' 
      else
        gridtype = 'REG'
      end if

      pi = 4.*atan(1.)
      hpi = pi/2.
      dtr = pi/180.
      aa = 6.37122e6
      omega = 7.292e-5
      fzero = omega+omega
      uzero = 40.
      pzero = 2.94e4
      alphad = 60.
      alpha = dtr*alphad
      indxn = (/((n,n=m,ntrunc),m=0,ntrunc)/)
      lap(:) = -real(indxn(:))*real(indxn(:)+1)/aa**2
c
      dt = 300.
      itmax = nint(86400.*5./dt)
      mprint = itmax/10

c
c     compute the derivative of the unrotated geopotential
c             p as a function of latitude
c
      nl = 91
      nlm1 = nl-1
      nlm2 = nl-2
      cfn = 1./nlm1
      dlath = pi/nlm1
      do 10 i=1,nlm2
      theta = i*dlath
      sth = sin(theta)
      cth = cos(theta)
      uhat = ui(uzero,hpi-theta)
      phlt(i) = cfn*cth*uhat*(uhat/sth+aa*fzero)
   10 continue
c
c     compute sine transform of the derivative of the geopotential
c     for the purpose of computing the geopotential by integration
c     see equation (3.9) in reference [1] above
c
      call sine(nlm2,phlt)
c
c     compute the cosine coefficients of the unrotated geopotential
c     by the formal integration of the sine series representation
c
      do 12 i=1,nlm2
      phlt(i) = -phlt(i)/i
   12 continue
c
c     phlt(i) contains the coefficients in the cosine series
c     representation of the unrotated geopotential that are used
c     below to compute the geopotential on the rotated grid.
c
c     compute the initial values of  east longitudinal
c     and latitudinal velocities u and v as well as the
c     geopotential p and coriolis f on the rotated grid.
c
      ca = cos(alpha)
      sa = sin(alpha)
      dlam = (pi+pi)/nlon
      dtheta = pi/(nlat-1)

      if (gridtype .eq. 'GAU') call gauinfo(gaulats,weights)

      do j=1,nlon
      lambda = (j-1)*dlam
      cl = cos(lambda)
      sl = sin(lambda)
      do i=1,nlat
c
c     lambda is longitude, theta is colatitude, and pi/2-theta is
c     latitude on the rotated grid. lhat and that are longitude
c     and colatitude on the unrotated grid. see text starting at
c     equation (3.10)
c
      if (gridtype .eq. 'GAU') then
         theta = hpi-gaulats(i)
      else
         theta = (i-1)*dtheta
      end if
      st = cos(theta)
      ct = sin(theta)
      sth = ca*st+sa*ct*cl
      cthclh = ca*ct*cl-sa*st
      cthslh = ct*sl
      lhat = atanxy(cthclh,cthslh)
      clh = cos(lhat)
      slh = sin(lhat)
      cth = clh*cthclh+slh*cthslh
      that = atanxy(sth,cth)
      uhat = ui(uzero,hpi-that)
      pxact(j,i) = cosine(that,nlm2,phlt)
      uxact(j,i) = uhat*(ca*sl*slh+cl*clh)
      vxact(j,i) = uhat*(ca*cl*slh*st-clh*sl*st+sa*slh*ct)
      f(j,i) = fzero*sth
      enddo
      enddo

      vmax = 0.
      pmax = 0.
      v2max = 0.
      p2max = 0.
      do j=1,nlat
      do i=1,nlon
         v2max = v2max+uxact(i,j)**2+vxact(i,j)**2
         p2max = p2max+pxact(i,j)**2
         vmax = amax1(abs(uxact(i,j)),abs(vxact(i,j)),vmax)
         pmax = amax1(abs(pxact(i,j)),pmax)
      enddo
      enddo
c
c     initialize first time step
c
      u = uxact
      v = vxact
      p = pxact
      ug = u
      vg = v
      pg = p

c  compute spectral coeffs of initial vrt,div,p.

      call getvrtdivspec(ug,vg,vrtnm,divnm,aa,gridtype)
      call grdtospec(pg,pnm,gridtype)

c==> time step loop.

      nnew = 1
      nnow = 2
      nold = 3
      do ncycle=0,itmax
 
      time = float(ncycle)*dt

C==> INVERSE TRANSFORM TO GET VORT AND PHIG ON GRID.
 
      call spectogrd(vrtnm,vrtg,gridtype)
      call spectogrd(pnm,pg,gridtype)
 
c==> compute u and v on grid from spectral coeffs. of vort and div.
 
      call getuv(vrtnm,divnm,ug,vg,aa,gridtype)

c==> compute error statistics.

      if (mod(ncycle,mprint) .eq. 0) then
      call spectogrd(divnm,divg,gridtype)
      u = ug
      v = vg
      p = pg
      htime = time/3600.
      write(*,390) ncycle,htime,dt,nlat,nlon,ntrunc,omega,pzero,
     1             uzero,alphad
  390 format(//' steady nonlinear rotated flow, test case 3',/
     1         ' cycle number              ',i10,
     2         ' model time in  hours      ',f10.2,/
     3         ' time step in seconds      ',f10.0,
     4         ' number of latitudes       ',i10,/
     5         ' number of longitudes      ',i10,
     6         ' max wave number           ',i10,/
     7         ' rotation rate        ',1pe15.6,
     8         ' mean height          ',1pe15.6,/
     9         ' maximum velocity     ',1pe15.6,
     1         ' tilt angle           ',1pe15.6)
      dvgm = 0.
      dvmax = 0.
      dpmax = 0.
      evmax = 0.0
      epmax = 0.0
      do j=1,nlat
      do i=1,nlon
         dvgm = amax1(dvgm,abs(divg(i,j)))
         dvmax = dvmax+(u(i,j)-uxact(i,j))**2+(v(i,j)-vxact(i,j))**2
         dpmax = dpmax+(p(i,j)-pxact(i,j))**2
         evmax = 
     *   amax1(evmax,abs(v(i,j)-vxact(i,j)),abs(u(i,j)-uxact(i,j)))
         epmax = amax1(epmax,abs(p(i,j)-pxact(i,j)))
      enddo
      enddo
      dvmax = sqrt(dvmax/v2max)
      dpmax = sqrt(dpmax/p2max)
      evmax = evmax/vmax
      epmax = epmax/pmax
      write(*,391) evmax,epmax,dvmax,dpmax,dvgm
  391 format(' max error in velocity',1pe15.6,
     +       ' max error in geopot. ',1pe15.6,/
     +       ' l2 error in velocity ',1pe15.6,
     +       ' l2 error in geopot.  ',1pe15.6,/
     +       ' maximum divergence   ',1pe15.6)
      end if

C==> COMPUTE RIGHT-HAND SIDES OF PROGNOSTIC EQNS. 
 
      scrg1(:,:) = ug(:,:)*(vrtg(:,:)+f(:,:))
      scrg2(:,:) = vg(:,:)*(vrtg(:,:)+f(:,:))
      call getvrtdivspec(scrg1,scrg2,ddivdtnm(:,nnew),
     *  dvrtdtnm(:,nnew),aa,gridtype)
      dvrtdtnm(:,nnew) = -dvrtdtnm(:,nnew)
      scrg1(:,:) = ug(:,:)*(pg(:,:)+pzero)
      scrg2(:,:) = vg(:,:)*(pg(:,:)+pzero)
      call getvrtdivspec(scrg1,scrg2,scrnm,dpdtnm(:,nnew),aa,gridtype)
      dpdtnm(:,nnew) = -dpdtnm(:,nnew)
      scrg1(:,:) = pg(:,:)+0.5*(ug(:,:)**2+vg(:,:)**2)
      call grdtospec(scrg1,scrnm,gridtype)
      ddivdtnm(:,nnew)=ddivdtnm(:,nnew)-lap(:)*scrnm(:)
 
c==> update vrt and div with third-order adams-bashforth.

c==> forward euler, then 2nd-order adams-bashforth time steps to start.

      if (ncycle .eq. 0) then
         dvrtdtnm(:,nnow) = dvrtdtnm(:,nnew)
         dvrtdtnm(:,nold) = dvrtdtnm(:,nnew)
         ddivdtnm(:,nnow) = ddivdtnm(:,nnew)
         ddivdtnm(:,nold) = ddivdtnm(:,nnew)
         dpdtnm(:,nnow) = dpdtnm(:,nnew)
         dpdtnm(:,nold) = dpdtnm(:,nnew)
      else if (ncycle .eq. 1) then
         dvrtdtnm(:,nold) = dvrtdtnm(:,nnew)
         ddivdtnm(:,nold) = ddivdtnm(:,nnew)
         dpdtnm(:,nold) = dpdtnm(:,nnew)
      end if
      vrtnm(:) = vrtnm(:) + dt*(
     *   (23./12.)*dvrtdtnm(:,nnew) - (16./12.)*dvrtdtnm(:,nnow)+
     *   (5./12.)*dvrtdtnm(:,nold) ) 
      divnm(:) = divnm(:) + dt*(
     *   (23./12.)*ddivdtnm(:,nnew) - (16./12.)*ddivdtnm(:,nnow)+
     *   (5./12.)*ddivdtnm(:,nold) ) 
      pnm(:) = pnm(:) + dt*(
     *   (23./12.)*dpdtnm(:,nnew) - (16./12.)*dpdtnm(:,nnow)+
     *   (5./12.)*dpdtnm(:,nold) ) 
 
C==> SWITCH INDICES.
 
      nsav1 = nnew
      nsav2 = nnow
      nnew = nold  
      nnow = nsav1
      nold = nsav2

      enddo

c garbage collection.

      call cleanup(gridtype)

      stop
      end
c
      function ui(amp,thetad)
c
c     computes the initial unrotated longitudinal velocity
c     see section 3.3.
c
      pi=4.*atan(1.)
      thetab=-pi/6.
      thetae= pi/2.
      xe=3.e-1
      x =xe*(thetad-thetab)/(thetae-thetab)
      ui = 0.
      if(x.le.0. .or. x.ge.xe) return
      ui=amp*exp(-1./x-1./(xe-x)+4./xe)
      return
      end
c
      function atanxy(x,y)
      atanxy = 0.
      if(x.eq.0. .and. y.eq.0.) return
      atanxy = atan2(y,x)
      return
      end
      subroutine sine(n,x)
c
c     computes the sine transform
c
      dimension x(n),w(n)
      arg = 4.*atan(1.)/(n+1)
      do 10 j=1,n
      w(j) = 0.
      do 10 i=1,n
      w(j) = w(j)+x(i)*sin(i*j*arg)
   10 continue
      do 15 i=1,n
      x(i) = 2.*w(i)
   15 continue
      return
      end
c
      function cosine(theta,n,cf)
c
c     computes the cosine transform
c
      dimension cf(n)
      cosine = 0.
      do 10 i=1,n
      cosine = cosine+cf(i)*cos(i*theta)
   10 continue
      return
      end
