module baro_integrals

     use precisions, only: wp, cwp
     use dispmodule
     use my_sparse_aggregate
     use my_sparse
     use my_trigs
     use save_load!, only:load_alloc_vector, save_vector
     use sal
     use my_blas

     type fric_breakdown
          real(wp) :: shallow, deep, coastal, open, total
     end type
     type domain_integrals
          real(wp) :: mass, ke, pe, d, dsal, df
          type(fric_breakdown) :: dit, dit_m, dit_p, dbl
     end type

contains

!**********************************************************************************************************
!**********************************************************************************************************
subroutine calc_convergence(nu, nv, nh, itn, P, cpts, delta, dir_cols, dir_sols)

    implicit none

     type(params), intent(in)	:: P

     integer, intent(in) 		:: nu, nv, nh, itn
     character(len=*), intent(in)	:: cpts, dir_cols, dir_sols
     real(wp), allocatable     	:: delta(:, :, :) !delta(itn,type,ccpt)
     character(len=*), parameter:: word = 'Convg'
     character(len=len(word)), allocatable 	:: words(:)
     character(len=1), allocatable 	:: tabs(:)

	 character(len=2)		:: cpt
     complex(cwp)			:: u(nu), v(nv), h(nh)
     complex(cwp)			:: u_new(nu), v_new(nv), h_new(nh)
     real(wp), allocatable	:: ta_u(:),ta_v(:),ta_h(:), dA_u(:),dA_v(:),dA_h(:), dh(:)
     integer 				:: istat, ccpt, ncpts, coor
     type(disp_settings) 	:: ds


	coor = P%coor
	ncpts=len(cpts)/2
	allocate(words(ncpts+2),tabs(ncpts+2), stat = istat)
	tabs = char(9)
	words(1) = word
	words(2) = '-----'
	do ccpt = 1, ncpts
		cpt=cpts(2*ccpt-1:2*ccpt)
		words(ccpt+2) = cpt
	enddo

	  do ccpt = 1, ncpts

	     cpt=cpts(2*ccpt-1:2*ccpt)

	!    %==========================
	!    % calculate maximum errors
	!    %==========================
	    if ((itn == 1).and.(P%load_sol==0)) then
	          u=0.
	          v=0.
	          h=0.
	    else
	          call load_vector(u, dir_sols // cpt // '_u' // '_m0_old' // '.dat')
	          call load_vector(v, dir_sols // cpt // '_v' // '_m0_old' // '.dat')
	          call load_vector(h, dir_sols // cpt // '_h' // '_m0_old' // '.dat')
	    endif

	          call load_vector(u_new, dir_sols // cpt // '_u' // '_m0' // '.dat')
	          call load_vector(v_new, dir_sols // cpt // '_v' // '_m0' // '.dat')
	          call load_vector(h_new, dir_sols // cpt // '_h' // '_m0' // '.dat')

		allocate(dA_u(nu),dA_v(nv),dA_h(nh), dh(nh), stat = istat)
		call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
		call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')
		call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')
			if (coor == 1) then ! Mercator
				dA_u= 1/cosh(ta_u)**2
				dA_v= 1/cosh(ta_v)**2
				dA_h= 1/cosh(ta_h)**2
			elseif (coor == 2) then ! lat-lon
				dA_u= cos(ta_u)
				dA_v= cos(ta_v)
				dA_h= cos(ta_h)
			endif
		deallocate(ta_h)

	    delta(itn,1,ccpt) = maxval(abs(u_new-u))
	    delta(itn,2,ccpt) = maxval(abs(v_new-v))
	    dh = abs(h_new-h)
	    delta(itn,3,ccpt) = maxval(dh)
	    delta(itn,4,ccpt) = vec_vec_dot(dA_h,dh)/sum(dA_h)
	    dh=abs( abs(h_new) - abs(h) )
	    delta(itn,5,ccpt)=maxval(dh)
	    delta(itn,6,ccpt)=vec_vec_dot(dA_h,dh)/sum(dA_h)

	    delta(itn,7,ccpt) = vec_vec_dot(abs( (h_new-h)/(h_new*log2int(abs(h_new)>1e-6) &
	    												+ 1.*log2int(abs(h_new)<=1e-6)) ), dA_h)/sum(dA_h)
	    delta(itn,8,ccpt) = .5*(  vec_vec_dot(abs( (u_new-u)/(u_new*log2int(abs(u_new)>1e-6) &
	    												+ 1.*log2int(abs(u_new)<=1e-6)) ), dA_u)/sum(dA_u) + &
	    						  vec_vec_dot(abs( (v_new-v)/(v_new*log2int(abs(v_new)>1e-6) &
	    												+ 1.*log2int(abs(v_new)<=1e-6)) ), dA_v)/sum(dA_v)  )

	  enddo
	    call save_matrix_3d(delta, dir_sols // 'delta.dat')

	call disp(char(9)//'-----------------------------------------------------------')

	ds = disp_get()
	  CALL DISP(tabs, advance='no')
	  CALL DISP_SET(STYLE='underline', SEP=char(9), MATSEP=' | ')!, ORIENT = 'ROW')

	  CALL DISP(words, TRIM='yes',advance='no')
	  CALL DISP('du', delta(itn,1,:), FMT='F7.2',TRIM='yes', advance='no')
	  CALL DISP('dv', delta(itn,2,:), FMT='F7.2',TRIM='yes', advance='no')
	  CALL DISP('duv%', 100*delta(itn,8,:), FMT='F6.1',TRIM='yes', advance='no')
	  if (P%cvg_scheme == 1) then
		  CALL DISP('dh', delta(itn,3,:), FMT='F7.2',TRIM='yes', advance='no')
		  CALL DISP('dhbar', delta(itn,4,:), FMT='F6.3',TRIM='yes', advance='no')
	  else
		  CALL DISP('dh', delta(itn,5,:), FMT='F6.3',TRIM='yes', advance='no')
		  CALL DISP('dhbar', delta(itn,6,:), FMT='F6.3',TRIM='yes', advance='no')
	  endif
	  CALL DISP('dh%', 100*delta(itn,7,:), FMT='F6.1',TRIM='yes', advance='no')

	  CALL DISP(char(9))
!	  CALL DISP('du'//char(9)//'dv'//char(9)//'dh'//char(9)//'dhbar'//char(9)//'dh'//char(9)//'dhbar',&
!	  			transpose(delta(itn,:,:)), FMT='F7.3',TRIM='yes')
	call disp_set(ds)
	call disp(char(9)//'-----------------------------------------------------------')

end subroutine calc_convergence
!**********************************************************************************************************

!**********************************************************************************************************
subroutine show_domain_integrals(cpt, di, di_log)

implicit none

     type(domain_integrals), intent(in) :: di
     character(len=2), intent(in)		:: cpt
     integer, optional, intent(in)		:: di_log

     type(disp_settings) ds

		ds = disp_get()
		call tostring_set(rfmt='F12.1')

		call disp('-----------------------------------------------------------------------------')
	    call disp(cpt//':'//char(9)//'KE = '//tostring(di%ke/1d15, 'F6.1')//' PJ, PE = '//tostring(di%pe/1d15, 'F6.1')//&
	    			' PJ, D = '//tostring(di%d/1d12, 'F6.3')//&
	              ' TW. D_SAL = '//tostring(di%dsal/1d9, 'F6.3')//' GW, D_f = '//tostring(di%df/1d9, 'F6.3')//' GW.')

	    call disp(char(9)//'D_BL = '//tostring(di%dbl%total/1d12, 'F6.3')//' TW, D_IT = '//tostring(di%dit%total/1d12, 'F6.3')//&
	    		  ' TW: D_ITP = '//tostring(di%dit_p%total/1d12, 'F6.3')//' TW, D_ITM = '//tostring(di%dit_m%total/1d9, 'F6.1')//' GW.' )
		call disp(char(9), '---------------------------------------------------------------------')

		  CALL DISP([char(9),char(9),char(9),char(9),char(9),char(9)], advance='no')
		  CALL DISP_SET(STYLE='underline', SEP=char(9), MATSEP=' | ')

		  CALL DISP('****', ['D_BL','DITP','DITM', 'D   '], advance='no')
		  CALL DISP('Total', [di%dbl%total/1d12,di%dit_p%total/1d12,di%dit_m%total/1d12,(di%dbl%total + di%dit%total)/1d12], FMT='F5.3',&
		  																					TRIM='yes', advance='no')
		  CALL DISP('Shallow', [di%dbl%shallow/1d12,di%dit_p%shallow/1d12,di%dit_m%shallow/1d12,(di%dbl%shallow + di%dit%shallow)/1d12], &
		  																		FMT='F5.3', TRIM='yes', advance='no')
		  CALL DISP('Deep', [di%dbl%deep/1d12,di%dit_p%deep/1d12,di%dit_m%deep/1d12,(di%dbl%deep + di%dit%deep)/1d12], &
		  																		FMT='F5.3', TRIM='yes', advance='no')
		  if ( abs((di%dbl%coastal + di%dbl%open)/di%dbl%total - 1) < 0.01 ) then ! check that they've been calculated
		  CALL DISP('Coastal', [di%dbl%coastal/1d12,di%dit_p%coastal/1d12,di%dit_m%coastal/1d12,(di%dbl%coastal + di%dit%coastal)/1d12], &
		  																		FMT='F5.3', TRIM='yes', advance='no')
		  CALL DISP('Open', [di%dbl%open/1d12,di%dit_p%open/1d12,di%dit_m%open/1d12,(di%dbl%open + di%dit%open)/1d12], FMT='F5.3', &
		  																					TRIM='yes', advance='no')
		  endif
		  CALL DISP(char(9))
!		  CALL DISP('Total'//char(9)//'Shallow'//char(9)//'Deep'//char(9)//'Coastal'//char(9)//'Open', &
!			transpose(reshape([di%dbl%total/1d12, di%dbl%shallow/1d12, di%dbl%deep/1d12, di%dbl%coastal/1d12, di%dbl%open/1d12, &
!		  			  di%dit%total/1d12, di%dit%shallow/1d12, di%dit%deep/1d12, di%dit%coastal/1d12, di%dit%open/1d12, &
!		  			  (di%dbl%total + di%dit%total)/1d12, (di%dbl%shallow + di%dit%shallow)/1d12, (di%dbl%deep + di%dit%deep)/1d12, &
!					  (di%dbl%coastal + di%dit%coastal)/1d12,(di%dbl%open + di%dit%open)/1d12], (/ 5, 3 /))), FMT='F5.3',TRIM='yes')

		call disp_set(ds)

!	  CALL DISP(words, TRIM='yes',advance='no')


end subroutine show_domain_integrals
!**********************************************************************************************************

subroutine log_di(cpt, itn, di, di_name, dir)

implicit none

	type(domain_integrals), intent(in) :: di
	character(len = *) :: dir, di_name!, cpts
	integer :: itn, c, fh = 16
	character(len=2)	:: cpt

    type(disp_settings) ds


		!*********************************
		! Save domain integrals in a log
		!*********************************
		open(fh, file=dir//di_name, ACCESS = 'APPEND')

	  call disp('==============',  UNIT=fh)
	  call disp(' Iteration ', itn,  UNIT=fh)
	  call disp('==============',  UNIT=fh)

		ds = disp_get()
		call tostring_set(rfmt='F12.1')

!do c = 1, len(cpts)/2
!		cpt = cpts(2*c-1:2*c)

	    call disp(cpt//':'//char(9)//'KE = '//tostring(di%ke/1d15, 'F6.1')//' PJ, PE = '//tostring(di%pe/1d15, 'F6.1')//&
	    			' PJ, D = '//tostring(di%d/1d12, 'F6.3')//' TW.'//&
	              ' D_SAL = '//tostring(di%dsal/1d9, 'F6.3')//' GW, D_f = '//tostring(di%df/1d9, 'F6.3')//' GW.', UNIT=fh)

	    call disp(char(9)//'D_BL = '//tostring(di%dbl%total/1d12, 'F6.3')//' TW, D_IT = '//tostring(di%dit%total/1d12, 'F6.3')//&
	' TW: D_ITP = '//tostring(di%dit_p%total/1d12, 'F6.3')//' TW, D_ITM = '//tostring(di%dit_m%total/1d9, 'F6.1')//' GW.', UNIT=fh)
		call disp(char(9), '---------------------------------------------------------------------', UNIT=fh)

		  CALL DISP([char(9),char(9),char(9),char(9),char(9)], advance='no', UNIT=fh)
		  CALL DISP_SET(STYLE='underline', SEP=char(9), MATSEP=' | ', UNIT=fh)

		  CALL DISP('****', ['D_BL','DITP','DITM', 'D   '], advance='no', UNIT=fh)
		  CALL DISP('Total', [di%dbl%total/1d12,di%dit_p%total/1d12,di%dit_m%total/1d12,(di%dbl%total + di%dit%total)/1d12], FMT='F5.3',&
		  																					TRIM='yes', advance='no', UNIT=fh)
		  CALL DISP('Shallow', [di%dbl%shallow/1d12,di%dit_p%shallow/1d12,di%dit_m%shallow/1d12,(di%dbl%shallow + di%dit%shallow)/1d12], &
		  																		FMT='F5.3', TRIM='yes', advance='no', UNIT=fh)
		  CALL DISP('Deep', [di%dbl%deep/1d12,di%dit_p%deep/1d12,di%dit_m%deep/1d12,(di%dbl%deep + di%dit%deep)/1d12], &
		  																		FMT='F5.3', TRIM='yes', advance='no', UNIT=fh)
		  if ( abs((di%dbl%coastal + di%dbl%open)/di%dbl%total - 1) < 0.01 ) then ! check that they've been calculated
		  CALL DISP('Coastal', [di%dbl%coastal/1d12,di%dit_p%coastal/1d12,di%dit_m%coastal/1d12,(di%dbl%coastal + di%dit%coastal)/1d12], &
		  																		FMT='F5.3', TRIM='yes', advance='no', UNIT=fh)
		  CALL DISP('Open', [di%dbl%open/1d12,di%dit_p%open/1d12,di%dit_m%open/1d12,(di%dbl%open + di%dit%open)/1d12], FMT='F5.3', &
		  																					TRIM='yes', advance='no', UNIT=fh)
		  endif
		  CALL DISP(char(9), UNIT=fh)
!enddo

		call disp_set(ds)

		  close(fh)

end subroutine log_di

!**********************************************************************************************************
subroutine save_di(cpts, di, di_name, dir_sols)

implicit none

	type(domain_integrals) :: di(:)
	character(len = *) :: dir_sols, di_name, cpts
	integer :: c, fh = 15

!*********************************
! Save domain integrals in a file
!*********************************
  open(fh, file=dir_sols//di_name, action='write', status='replace')

  do c = 1, len(cpts)/2
  	write(fh, '(a)', advance='no') cpts(2*c-1:2*c)
  	write(fh, *) di(c)
  enddo

  close(fh)

end subroutine save_di


subroutine baro_domain_integrals(cpt, P, nu, nv, nh, itm_active, dir_grid, dir_cols, dir_mats, dir_sols, di)

!% calculates the (complex) equilibrium tide heq
!% on the h grid, for the specified domain
implicit none

     type(domain_integrals) :: di
     type(params)	:: P

     integer, intent(in)  :: nu, nv, nh, itm_active

     real(wp)  :: dA_u(nu), dA_v(nv), dA_h(nh)

     real(wp)	:: g, re, rhoe, rhoo, latP, lonP, cdg, Q, ubar, beta0
     integer	:: fr_scheme, sal_scheme, itd_scheme, cooruv

!     integer         :: hp(nh, 4), up(nu, 4), vp(nv, 4)
     real(wp)        :: ta_u(nu), ta_v(nv)!, ta_h(nh)
     real(wp)        :: H_h(nh), H_u(nu), H_v(nv)
     complex(cwp)    :: u(nu), v(nv), h(nh)!, u_l(nu), v_l(nv), h_l(nh)


     real(wp)         :: temp1(nu), temp2(nv), temp3(nh)
     complex(cwp)     :: Du(nu), Dv(nv)!, Du_l(nu), Dv_l(nv)
     complex(cwp), allocatable, dimension(:) :: tmp_u, tmp_v, tmp_u1, tmp_v1
     integer, allocatable, dimension(:) :: ufactor,vfactor
	 integer, allocatable	:: I_u(:), I_v(:)

!     character(len=*), intent(in)   :: cpts
     character(len=2)               :: cpt

	 logical			:: file_exist_u, file_exist_v
     integer            :: istatus

     type (triplet) :: v2uf, u2vf
     type (triplet) :: DUU, DVU, DUV, DVV
     complex(cwp)   :: hsal(nh)
     complex(cwp), allocatable, dimension(:) :: heq
     type (tide_params) :: tidal_pars

     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols
!==========================================================================================
! SHORTCUTS
!==========================================================================================
latP = P%latP
lonP = P%lonP
g = P%g
re = P%re
rhoe = P%rhoe
rhoo = P%rhoo
cdg = P%cd
ubar = P%ubar
Q = P%Qbar
beta0 = P%beta0

fr_scheme = P%fr_scheme
sal_scheme = P%sal_scheme
itd_scheme = P%itd_scheme
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
!%===============
!% load some grids
!%===============
!     call load_mat_int(up, dir_cols // 'up.dat')
!     call load_mat_int(vp, dir_cols // 'vp.dat')
!     call load_mat_int(hp, dir_cols // 'hp.dat')
!     call load_vector(ta_h, dir_cols // 'ta_h.dat')
     call load_vector(ta_u, dir_cols // 'ta_u.dat')
     call load_vector(ta_v, dir_cols // 'ta_v.dat')

     call load_vector(H_h, dir_cols // 'H_h.dat')
     call load_vector(H_u, dir_cols // 'H_u.dat')
     call load_vector(H_v, dir_cols // 'H_v.dat')
!%==============================
!% load solutions for every cpt
!%==============================
     call load_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
     call load_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
     call load_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

!%==============================
!% metric vectors for integration
!%==============================
    call get_dA(dA_u, dA_v, dA_h, re, P%coor, dir_grid, dir_cols)

!%==================================
!% calculate basic Ocean statistics
!%==================================
di%mass = rhoo*vec_vec_dot(dA_h, H_h)

	!***********************************************************************************
	! KINETIC ENERGY
	!***********************************************************************************
	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			temp1=0.25*rhoo*(abs(u)**2 * cosh(ta_u)**2)/H_u
			temp2=0.25*rhoo*(abs(v)**2 * cosh(ta_v)**2)/H_v
		case (22) ! LAT-LON, TRANSPORT
			temp1=0.25*rhoo*(abs(u)**2)/H_u
			temp2=0.25*rhoo*(abs(v)**2)/H_v
	end select

  di%ke = vec_vec_dot(dA_u, temp1) + vec_vec_dot(dA_v, temp2)

	!***********************************************************************************
	! POTENTIAL ENERGY
	!***********************************************************************************
	temp3 = 0.25*g*rhoo*abs(h)**2
	di%pe = vec_vec_dot(dA_h, temp3)

	!***********************************************************************************
	! calculate the globally integrated total dissipation (W)
	!***********************************************************************************
     call calc_heq_h(heq, cpt, nh, latP, lonP, P%coor, dir_cols, dir_grid)
     tidal_pars = get_pars(cpt)

     temp3 = 0.5*g*rhoo*tidal_pars%omega0*real( i*heq*conjg(h), kind=wp )
     di%d = vec_vec_dot(dA_h, temp3)
     deallocate(heq)

	!***********************************************************************************
	! calculate the globally integrated SAL dissipation (W)
	!***********************************************************************************
	if ((sal_scheme == 0).or.(itd_scheme==20)) then
	    hsal = 0
	elseif (sal_scheme == 1) then
	    hsal = beta0*h
	else
		call load_vector(hsal, dir_sols //'temp/'// cpt // '_hsal' // '.dat')
	endif

     temp3 = -0.5*g*rhoo*tidal_pars%omega0*real( i*hsal*conjg(h), kind=wp )
     di%dsal = vec_vec_dot(dA_h, temp3)

	!***********************************************************************************
	! now analyse the frictional dissipation (W/m^2):
	!***********************************************************************************
	if (fr_scheme == 0) then
		Du = 0
		Dv = 0
	elseif (fr_scheme == 1) then
		select case (cooruv)
			case (12) ! MERCATOR, TRANSPORT
				Du=cdg*ubar*cosh(ta_u)*u/H_u
				Dv=cdg*ubar*cosh(ta_v)*v/H_v
		    case (22) ! LAT-LON, TRANSPORT
				Du=cdg*ubar*u/H_u
				Dv=cdg*ubar*v/H_v
		end select
	elseif (fr_scheme == 2) then
		select case (cooruv)
			case (12) ! MERCATOR, TRANSPORT
		        Du=cdg*Q*cosh(ta_u)*u/(H_u**2)
		        Dv=cdg*Q*cosh(ta_v)*v/(H_v**2)
		    case (22) ! LAT-LON, TRANSPORT
		        Du=cdg*Q*u/(H_u**2)
		        Dv=cdg*Q*v/(H_v**2)
		end select
!	some kind of test
!     call load_vector(u_l, dir_sols // cpt // '_u' // '_m0' // '.dat')
!     call load_vector(v_l, dir_sols // cpt // '_v' // '_m0' // '.dat')
!          Du_l=cdg*Q*cosh(ta_u)*u_l/(H_u**2)
!          Dv_l=cdg*Q*cosh(ta_v)*v_l/(H_v**2)
!     temp1=0.5*rhoo*real(conjg(u_l)*Du_l, kind=wp)*cosh(ta_u)/H_u
!     temp2=0.5*rhoo*real(conjg(v_l)*Dv_l, kind=wp)*cosh(ta_v)/H_v
	elseif (fr_scheme == 3) then
		call load_vector(Du, dir_sols // cpt // '_Du' // '_m0_drag' // '.dat')
		call load_vector(Dv, dir_sols // cpt // '_Dv' // '_m0_drag' // '.dat')
	endif

	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			temp1=0.5*rhoo*real(conjg(u)*Du, kind=wp)*cosh(ta_u)/H_u
			temp2=0.5*rhoo*real(conjg(v)*Dv, kind=wp)*cosh(ta_v)/H_v
		case (22) ! LAT-LON, TRANSPORT
			temp1=0.5*rhoo*real(conjg(u)*Du, kind=wp)/H_u
			temp2=0.5*rhoo*real(conjg(v)*Dv, kind=wp)/H_v
	end select

 di%dbl%total = vec_vec_dot(dA_u, temp1) + vec_vec_dot(dA_v, temp2)
! using count(reshape(..., [n, 1]), 2) to convert logical to integer
 di%dbl%shallow = vec_vec_dot( dA_u, temp1*count(reshape(H_u <= P%sh_depth, [nu, 1]), 2) ) + &
 				  vec_vec_dot( dA_v, temp2*count(reshape(H_v <= P%sh_depth, [nv, 1]), 2) )
 di%dbl%deep = vec_vec_dot( dA_u, temp1*count(reshape(H_u > P%sh_depth, [nu, 1]), 2) ) + &
 				  vec_vec_dot( dA_v, temp2*count(reshape(H_v > P%sh_depth, [nv, 1]), 2) )

!	Calculate the coordinates of the coastal areas
if (P%coast_flag == 1) then
	inquire( file=dir_cols//'temp/'//'coast_ind_u.dat', exist=file_exist_u )
	inquire( file=dir_cols//'temp/'//'coast_ind_v.dat', exist=file_exist_v )
	if ( (.not. file_exist_u) .or. (.not. file_exist_v) ) then
		call calc_open_ocean(I_u,I_v, real(P%coast_dist, wp), P%coor, re, nu, nv, nh, dir_grid, dir_cols, dir_mats, P%lib)
!	Save to file if not kept in memory
		call save_vector(I_u, dir_cols//'temp/'//'coast_ind_u.dat')
		call save_vector(I_v, dir_cols//'temp/'//'coast_ind_v.dat')
	else
!	Load from file if not kept in memory
		call load_alloc_vector(I_u, dir_cols//'temp/'//'coast_ind_u.dat')
		call load_alloc_vector(I_v, dir_cols//'temp/'//'coast_ind_v.dat')
	endif

	 di%dbl%coastal = vec_vec_dot(dA_u, temp1*I_u) + vec_vec_dot(dA_v, temp2*I_v)
	 di%dbl%open = vec_vec_dot(dA_u, temp1*(1-I_u)) + vec_vec_dot(dA_v, temp2*(1-I_v))
else
	 di%dbl%coastal = 0._dp
	 di%dbl%open = 0._dp
endif

!  %====================================================
!  % now analyse the internal tide dissipation (W/m^2):
!  %====================================================
!	1. Parametrized ITD
    if ( (itd_scheme == 1).or.((itd_scheme >= 2).and.(itm_active==0)) &
			.or. ((itd_scheme == 3).and.(itm_active==1)) ) then

		call load_alloc_sparse(DUU, dir_mats // 'DUU.dat')
		call load_alloc_sparse(DUV, dir_mats // 'DUV.dat')
		call load_alloc_sparse(DVU, dir_mats // 'DVU.dat')
		call load_alloc_sparse(DVV, dir_mats // 'DVV.dat')
		call load_alloc_vector(ufactor, dir_cols // 'ufactor_'//cpt//'.dat')
		call load_alloc_vector(vfactor, dir_cols // 'vfactor_'//cpt//'.dat')

		allocate(tmp_u(nu), tmp_v(nv), stat = istatus)
		allocate(tmp_u1(nu), tmp_v1(nv), stat = istatus)

		call coo_vec_mul(DUU, u, P%lib, tmp_u)
		call coo_vec_mul(DVU, v, P%lib, tmp_u1)
		call coo_vec_mul(DUV, u, P%lib, tmp_v)
		call coo_vec_mul(DVV, v, P%lib, tmp_v1)

        Du = ufactor * (tmp_u + tmp_u1)/tidal_pars%omega0
        Dv = vfactor * (tmp_v + tmp_v1)/tidal_pars%omega0

        deallocate(tmp_u, tmp_v, tmp_u1, tmp_v1)
    else
		Du = 0.
		Dv = 0.
	end if

	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			temp1=0.5*rhoo*real(conjg(u)*Du, kind=wp)*cosh(ta_u)/H_u
			temp2=0.5*rhoo*real(conjg(v)*Dv, kind=wp)*cosh(ta_v)/H_v
		case (22) ! LAT-LON, TRANSPORT
			temp1=0.5*rhoo*real(conjg(u)*Du, kind=wp)/H_u
			temp2=0.5*rhoo*real(conjg(v)*Dv, kind=wp)/H_v
	end select

	di%dit_p%total = vec_vec_dot(dA_u, temp1) + vec_vec_dot(dA_v, temp2)
! using count(reshape(..., [n, 1]), 2) to convert logical to integer
 di%dit_p%shallow = vec_vec_dot( dA_u, temp1*count(reshape(H_u <= P%sh_depth, [nu, 1]), 2) ) + &
 				  vec_vec_dot( dA_v, temp2*count(reshape(H_v <= P%sh_depth, [nv, 1]), 2) )
 di%dit_p%deep = vec_vec_dot( dA_u, temp1*count(reshape(H_u > P%sh_depth, [nu, 1]), 2) ) + &
 				  vec_vec_dot( dA_v, temp2*count(reshape(H_v > P%sh_depth, [nv, 1]), 2) )
if (P%coast_flag == 1) then
	 di%dit_p%coastal = vec_vec_dot(dA_u, temp1*I_u) + vec_vec_dot(dA_v, temp2*I_v)
	 di%dit_p%open = vec_vec_dot(dA_u, temp1*(1-I_u)) + vec_vec_dot(dA_v, temp2*(1-I_v))
endif

!	do another calculation for modelled-only ITD
	if ((itd_scheme >= 2).and.(itm_active==1))  then
		call load_vector(Du, dir_sols//'itm/' // cpt // '_Du' // '_m0_itd_drag'// '.dat')
		call load_vector(Dv, dir_sols//'itm/' // cpt // '_Dv' // '_m0_itd_drag'// '.dat')
    else
		Du = 0.
		Dv = 0.
	end if

	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			temp1=0.5*rhoo*real(conjg(u)*Du, kind=wp)*cosh(ta_u)/H_u
			temp2=0.5*rhoo*real(conjg(v)*Dv, kind=wp)*cosh(ta_v)/H_v
		case (22) ! LAT-LON, TRANSPORT
			temp1=0.5*rhoo*real(conjg(u)*Du, kind=wp)/H_u
			temp2=0.5*rhoo*real(conjg(v)*Dv, kind=wp)/H_v
	end select

 di%dit_m%total = vec_vec_dot(dA_u, temp1) + vec_vec_dot(dA_v, temp2)
! using count(reshape(..., [n, 1]), 2) to convert logical to integer
 di%dit_m%shallow = vec_vec_dot( dA_u, temp1*count(reshape(H_u <= P%sh_depth, [nu, 1]), 2) ) + &
 				  vec_vec_dot( dA_v, temp2*count(reshape(H_v <= P%sh_depth, [nv, 1]), 2) )
 di%dit_m%deep = vec_vec_dot( dA_u, temp1*count(reshape(H_u > P%sh_depth, [nu, 1]), 2) ) + &
 				  vec_vec_dot( dA_v, temp2*count(reshape(H_v > P%sh_depth, [nv, 1]), 2) )
if (P%coast_flag == 1) then
	 di%dit_m%coastal = vec_vec_dot(dA_u, temp1*I_u) + vec_vec_dot(dA_v, temp2*I_v)
	 di%dit_m%open = vec_vec_dot(dA_u, temp1*(1-I_u)) + vec_vec_dot(dA_v, temp2*(1-I_v))
endif

! sum modeled and parametrized drag
di%dit%total = di%dit_m%total + di%dit_p%total
di%dit%shallow = di%dit_m%shallow + di%dit_p%shallow
di%dit%deep = di%dit_m%deep +  di%dit_p%deep
if (P%coast_flag == 1) then
	 di%dit%coastal = di%dit_m%coastal + di%dit_p%coastal
	 di%dit%open = di%dit_m%open + di%dit_p%open
endif


!  %===========================================================
!  % finally check for Coriolis energy input (should be zero)
!  %===========================================================
	call load_alloc_sparse(v2uf, dir_mats // 'v2uf.dat')
	call load_alloc_sparse(u2vf, dir_mats // 'u2vf.dat')

	allocate(tmp_u(nu), tmp_v(nv), stat = istatus)
	call coo_vec_mul(v2uf, v, P%lib, tmp_u)
	call coo_vec_mul(u2vf, u, P%lib, tmp_v)

	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			temp1 = -0.5*rhoo*real(conjg(u)*tmp_u, kind=wp)*cosh(ta_u)**2/H_u
			temp2 = 0.5*rhoo*real(conjg(v)*tmp_v, kind=wp)*cosh(ta_v)**2/H_v
		case (22) ! LAT-LON, TRANSPORT
			temp1 = -0.5*rhoo*real(conjg(u)*tmp_u, kind=wp)/H_u
			temp2 = 0.5*rhoo*real(conjg(v)*tmp_v, kind=wp)/H_v
	end select

deallocate(tmp_u, tmp_v)
call dealloc_sparse(u2vf)
call dealloc_sparse(v2uf)

    di%df = vec_vec_dot(dA_u, temp1) + vec_vec_dot(dA_v, temp2)
 !call disp('temp1: ', temp1, ADVANCE='NO')
 !call disp('temp1: ', temp2)

end subroutine baro_domain_integrals

!**********************************************************************************************************
!**********************************************************************************************************

subroutine get_dA(dA_u, dA_v, dA_h, re, coor, dir_grid, dir_cols)

implicit none

     real(wp)           :: dph, dta, dth

     real(wp), intent(in) :: re
     integer, intent(in)  :: coor

     real(wp) :: dA_u(:), dA_v(:), dA_h(:)
     real(wp), allocatable :: ta_h(:), ta_u(:), ta_v(:)

     real(wp), allocatable :: ph_vg(:), ta_ug(:), th_ug(:)

!     integer ::      istatus
     character(len = *) :: dir_grid, dir_cols


     call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
     call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')
     call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')

!%===================
!% calculate metrics
!%===================
	select case (coor)
	!***********************************************************************************
		case (1) ! MERCATOR
			call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
			call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
			dph=ph_vg(2)-ph_vg(1)
			dta=ta_ug(2)-ta_ug(1)
			deallocate(ta_ug, ph_vg)

			dA_u= (re**2)*dta*dph/cosh(ta_u)**2
			dA_v= (re**2)*dta*dph/cosh(ta_v)**2
			dA_h= (re**2)*dta*dph/cosh(ta_h)**2
		case (2) ! LAT-LON
			call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
			call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
			dph=ph_vg(2)-ph_vg(1)
			dth=th_ug(2)-th_ug(1)
			deallocate(th_ug, ph_vg)

			dA_u= (re**2)*dth*dph*cos(ta_u)
			dA_v= (re**2)*dth*dph*cos(ta_v)
			dA_h= (re**2)*dth*dph*cos(ta_h)
	end select

	deallocate(ta_u, ta_v, ta_h)

end subroutine get_dA

!**********************************************************************************************************

subroutine calc_open_ocean(I_u, I_v, L, coor, re, nu, nv, nh, dir_grid, dir_cols, dir_mats, lib)

implicit none

     integer, intent(in)	:: coor, nu, nv, nh
     real(wp), intent(in)   :: re, L

     integer				:: nph, nth, cth, cph
     integer				:: nlat, nlon, istat, j, k
	 real(wp)           	:: dph, dth
	 integer, allocatable	:: ivals(:), jvals(:)
     integer, allocatable	:: I_u(:), I_v(:), I_h(:), I_hg(:, :)
     integer, allocatable 	:: hp(:, :)
     real(wp), allocatable	:: ph_vg(:), th_ug(:), H_hg(:, :), H_tmp(:,:), circle(:,:), tmp(:)

	 type(triplet)		:: h2u, h2v
     character(len = *) :: dir_grid, dir_cols, dir_mats
     character			:: lib

!	Load
	call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat', nph)
	call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat', nth)
	call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
	call load_alloc_matrix(hp, dir_cols // 'hp.dat')

	allocate(I_hg(nph, nth), stat=istat)
	I_hg = 0


do cth = 1, nth

!	dph, dth -- linear dims of a cell aroun point ( ph_vg(cth), th_ug(cth) )
	if (coor == 1) then ! Mercator (~ square cells)
		dph = 2*pi*re*cos(th_ug(cth))/nph
		nlon = L/dph ! does rounding
		dth = dph
		nlat = nlon
	else
		dth = 2*pi*re/nph
		nlat = L/dth ! does rounding
		dph = 2*pi*re*cos(th_ug(cth))/nph
		nlon = L/dph ! does rounding
	endif

		if (allocated(circle)) deallocate(circle)
			allocate(circle(-nlon:nlon, -nlat:nlat), stat = istat)
		circle = 0.
		if (allocated(H_tmp)) deallocate(H_tmp)
			allocate(H_tmp(-nlon:nlon, -nlat:nlat), stat = istat)
		H_tmp = 1 ! the fill-in value (for points outside of the cirlce, must be >0)

		if (allocated(ivals)) deallocate(ivals)
			allocate(ivals(-nlon:nlon), stat = istat)
		if (allocated(jvals)) deallocate(jvals)
			allocate(jvals(-nlat:nlat), stat = istat)
!%================================================================
!  % calc number approx number of grid points at this latitude...
!  % assume lat-lon grid for now
!%================================================================
	jvals = [ ( cth + j, j=-nlat,nlat ) ]
	where (jvals<1) jvals = 1
	where (jvals>nth) jvals = nth

	do cph = 1, nph
    	ivals = [ ( modulo(cph+k-1,nph)+1, k=-nlon, nlon ) ]
		do k = -nlon, nlon
			do j = -nlat, nlat
	    		circle(k, j) = (k*dph)**2 + (j*dth)**2
	    		if (circle(k, j) <= L**2) then
!	    			if (ivals(k)*jvals(j)==0) print *, ivals(k), jvals(j)
	    			H_tmp(k, j) = H_hg(ivals(k), jvals(j))
	    		endif
	    	enddo
	    enddo

!	    call disp(transpose(H_tmp))
!	    stop

	I_hg(cph,cth) = count(H_tmp<=0)
	enddo
!	call disp(minval(I_hg(:,cth)), advance = 'no')
enddo

	where (I_hg > 0) I_hg = 1
!	call save_matrix(I_hg, dir_mats // 'I_hg.dat')

    allocate(I_u(nu), I_v(nv), I_h(nh), stat = istat)
    I_u = 0
    I_v = 0
    I_h = 0

    do j = 1, nh
		I_h(j) = I_hg(hp(j,1),hp(j,2))
	end do

! load sparse matrices
    call load_alloc_sparse(h2u, dir_mats // 'h2u.dat')
	call load_alloc_sparse(h2v, dir_mats // 'h2v.dat')

    allocate(tmp(nu), stat = istat)
	call coo_vec_mul(h2u, real(I_h, kind = wp), lib, tmp)
	I_u = tmp
	deallocate(tmp)

    allocate(tmp(nv), stat = istat)
	call coo_vec_mul(h2v, real(I_h, kind = wp), lib, tmp)
	I_v = tmp
	deallocate(tmp)

    call deallocate_sparse(h2u)
    call deallocate_sparse(h2v)

!deallocate(hp, ph_vg, th_ug)

end subroutine calc_open_ocean

!**********************************************************************************************************

end module baro_integrals
