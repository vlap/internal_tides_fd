program baro_v1

use precisions, only: wp, cwp, dp
use my_paths
use baro_integrals
use baro_solver_mod
use itm_solver_mod
use generate_matrices
use generate_grid
use my_trigs
use my_sparse
use save_load
use control
use dispmodule


use itd

     implicit none

     type(params)	:: P, P_load
     type(grid_dims):: GD

	 real(dp), allocatable	:: N(:), z(:), H_h(:)
	 real(dp)				:: Nbar(1)
	 integer				:: nz
!     integer            :: nlon, nlat, nz, clon, clat
!     real(wp), allocatable	:: lon(:), lat(:), z(:)
!     real(wp), allocatable	:: N_topo_3d(:,:,:)



!	integer	::  j, m
!	integer	:: nmodes, norder
!	real(wp) :: hmax, hmin
!	real(wp), allocatable	:: c2n(:,:), wmodes(:,:)

	real(dp), allocatable	:: c2n_interp(:,:), dc2n_interp(:,:) ,xinterp(:)


!**********************************************************************************************************
!**********************************************************************************************************
	 logical			:: dir_e=.false., identical_P=.false.
     integer			:: ccpt, ncpts
     character(len=2)	:: cpt
     character(len=10)	:: command

     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max
     type(domain_integrals), allocatable :: di(:)

!********************************************************************

     call system('clear')
!%================================
!% load parameters of the problem:
!%================================
    call control_file(P_file, P)
	!	Quick consistency check
	call control_consistency(P)

!%========================
!% set up file structure:
!%========================
     call make_save_dirs(P%isave, save_gridid, P%itd_scheme>=2, len(out_dir), out_dir, &
     															dir_cols, dir_grid, dir_mats, dir_sols, dir_global)
	 call system('yes | cp control_file.txt ' //  dir_global)

!%================================
!% Generate global grid & matrices
!% 	or use previously generated
!%================================
if (len(trim(P%gridid)) < len('0000_00_00_00_00')) then
	!%================================================
	!% 1) Loading/preparing topography file
	!% 2) Set up the global C-grid and mapping matrices
	!%================================================
		call generate_global_grid(etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols, P, GD)
		! Write the parameters of the grid into a specified file
		call write_GD(dir_grid//GD_file, GD)
	!%================================================
	!% 1) Write H, ta on the u/v/h-grids
	!% 2) Generate and write sparse matrices
	!%================================================
		call generate_global_matrices(dir_grid, dir_cols, dir_mats, P, GD)
else
	! nph, latP and omega should remain the same as in that run
	inquire( file=out_dir//trim(P%gridid)//'/.', exist=dir_e )
	if  (.not. dir_e ) then
!		Use MATLAB functions to export the data in M-files into binaries
		call system('xterm -e matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
		'"addpath(genpath(''' // matlab_dir // ''')); save_grid_binary(''' // trim(P%gridid) // ''', '''//trim(model)// '''); exit;" ')
!		print *, ' matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
!		'"addpath(genpath(''' // matlab_dir // ''')); save_grid_binary(''2012_06_13__16_54'')" '

!		if the matlab script failed
		inquire( file=out_dir//trim(P%gridid)//'/.', exist=dir_e )
		if  (.not. dir_e ) then
			write(*,'(a)') 'Directory /'//trim(P%gridid)//' for the previously generated grid files doesn''t exist.'
			stop
		end if
	else
		!	Load the parameters file
		call control_file(out_dir // trim(P%gridid) // '/' // P_file, P_load)
		identical_P = cmpr_P_grids(P, P_load)
		if (.not.( identical_P )) then
			print *, "*************************************************************************************************"
			print *, "The parameters for the previously generated and current grid do not match exactly."
			print *, "Generating the grid from scratch..."
			print *, "*************************************************************************************************"
		endif
	end if

		! NOTE: hard linking creates problems when old dirs are removed.
		if (P%grid_link==1) then
			command = 'ln -s '
		else
			command = 'rsync -a '
		endif

	if ((.not. dir_e ).or.(identical_P)) then

		if (P%grid_link==1) then
			write(*, '("Linking to previously generated grid: /", a)')trim(P%gridid)
		else
			write(*, '("Copying previously generated grid: /", a)')trim(P%gridid)
		endif

	!%================================================
	!% Link to previously generated files (either exported with matlab, or generated earlier)
	!%================================================

	     call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/'//'*.dat ' // dir_grid)
	     call system(command // out_dir//trim(P%gridid)//'/'//'global/cols/'//'*.dat ' // dir_cols)

	!	if (P%grid_link==1) then
			call system(command  // out_dir//trim(P%gridid)//'/'//'global/mats/'//'u2*.dat ' // dir_mats)
			call system(command  // out_dir//trim(P%gridid)//'/'//'global/mats/'//'v2*.dat ' // dir_mats)
			call system(command  // out_dir//trim(P%gridid)//'/'//'global/mats/'//'h2*.dat ' // dir_mats)
	!		call system('rm -f ' // dir_mats // 'mat_*.dat')
	!	else
	!		call system(command // ' --exclude="mat_*.dat" ' // out_dir//trim(P%gridid)//'/'//'global/mats/*.dat ' // dir_mats)
	!	endif

	!     call system('find ' // out_dir//trim(P%gridid)//'/'//&
	!		 'global/mats/ -type f -not -iname ''mat_*.dat'' -exec cp ''{}'' '''//&
	!		 dir_mats//''' '';''')

		! Read the parameters of the grid into a specified file
	     call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/'//'*.txt ' // dir_grid)
	     call read_GD(dir_grid//GD_file, GD)
	     ! Update P to match params of the uploaded grid
	     P%nph  = GD%nph
	     P%coor = GD%coor
	     P%lonP = GD%lonP
	     P%latP = GD%latP
	     	! check if after the adjustment the grids are the same
			if (.not.( identical_P )) then
				identical_P = cmpr_P_grids(P, P_load)
			endif

	else
		!%================================================
		!% 1) Loading/preparing topography file
		!% 2) Set up the global C-grid and mapping matrices
		!%================================================
			call generate_global_grid(etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols, P, GD)
			! Write the parameters of the grid into a specified file
			call write_GD(dir_grid//GD_file, GD)
		!%================================================
		!% 1) Write H, ta on the u/v/h-grids
		!% 2) Generate and write sparse matrices
		!%================================================
			call generate_global_matrices(dir_grid, dir_cols, dir_mats, P, GD)
	endif

end if
!%================================================
!% Attempt to load the initial guess baro solution
!%================================================
if (P%load_sol >= 1) then
	call prepare_init_guess(P, GD, out_dir,P_file, GD_file, tidal_dir, dir_grid, dir_cols, dir_mats, dir_sols)
endif

!%============================
!% IT modelling. Cn and sponge
!%============================
if (P%itd_scheme >= 2) then

  call CPU_Time(T1)
	call system_clock ( wall_t1, clock_rate, clock_max )


    print *, ""

! If the cn was found the same way on the previous grid, just copy the results

! commented out ad-hoc not to recalculate cn for smoothed baro runs!
	identical_P=.true. ! SHOULD BE REMOVED!

	if (identical_P) then
		inquire( file=out_dir//trim(P%gridid)//'/global/grid/'//'/itm/c2n_h.dat', exist=dir_e )
		identical_P = cmpr_P_cn(P, P_load)

		if ( (dir_e).and.(identical_P) ) then
			write(*, '("Linking to c_n on the previously generated grid: ")', advance='no')
     		call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/itm/'//'c2n_u.dat ' // dir_grid//'/itm/')
     		call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/itm/'//'c2n_v.dat ' // dir_grid//'/itm/')
     		call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/itm/'//'c2n_h.dat ' // dir_grid//'/itm/')
     		call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/itm/'//'c2n_h_orig.dat ' // dir_grid//'/itm/')
     		call system(command // out_dir//trim(P%gridid)//'/'//'global/grid/itm/'//'dc2n_h.dat ' // dir_grid//'/itm/')
		endif
	endif

! Calculate phase speeds c_m and its derivatives
	if ( .not.((dir_e).and.(identical_P)) ) then
	!	hmin = 300._dp
	!	hmax = 11000._dp
		write(*, '("Solve the eigenv problem to find c_n: ")', advance='no')

		if ( (P%N_hor_var == 0).or.(P%cn_scheme==1) ) then
			!%================================================
			!% 1) Calculate average N(z) from WOA09
			!% 2) Solve the eigenv problem to find c_m(H(x))
			!%================================================
		    call calc_N_z_avrg(P, N_data_dir, dir_grid)
			call load_alloc_vector(z, dir_grid // 'z.dat')
			call load_alloc_vector(N, dir_grid // 'N_avrg.dat')

			if ((P%cn_scheme==1) ) then
				nz = size(z)
				Nbar = 1/(z(nz)-z(1)) * sum( (N(2:nz)+N(1:nz-1))/2 * (z(2:nz)-z(1:nz-1)) ) ! trapz rule averaging
			    call eval_c2n_of_H(P, P%N_hmin, P%N_hmax, z, Nbar, dir_grid, dir_cols)
			else
				call eval_c2n_of_H(P, P%N_hmin, P%N_hmax, z, N, dir_grid, dir_cols)
			endif
		elseif (P%N_hor_var == 1) then
			!%================================================
			!% 1) Typical atlas grid, like WOA09, is 360x180,
			!% then treat each grid cell as a subdomain where
			!% 2) Solve the eigenv problem to find c_m(H(x))
			!%================================================
			write(*, '("horizontally non-uniform subdomains...")')
			call eval_c2n_of_H_and_N(P, GD, N_data_dir, P%N_hmin, P%N_hmax, dir_grid, dir_cols, dir_mats)
			! slower version that solved the eigenproblem on both u,v-grids instead of interpolating from h-grid
!			call eval_c2n_of_H_and_N_old(P, N_data_dir, P%N_hmin, P%N_hmax, dir_grid, dir_cols)
		else
			print *, "Choose appropriate cn_scheme: 1 or 2"
			stop
		endif
		!%===================================
	endif

  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall: '&
						//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')
    print *, ""

	!%================================================
	!% Generate sponge. Note: frequency dependent
	!%================================================
	call generate_sponge(trim(P%cpts), P, GD, dir_grid, dir_cols, dir_mats)

 endif

!%================================================
!% 1) Collect all the sparse matrices together
!% 2) Write mat and bcdiag
!%================================================
	call write_baro_mats(N_data_dir, dir_cols, dir_grid, dir_mats, trim(P%cpts), P, GD)
write(*, '("----------------------------------------------------------------------------------------")')
!%================================================
!% Do the same for IT modelling matrices
!%================================================
if (P%itd_scheme >= 2) then
	call write_itm_mats(trim(P%cpts), P, GD, dir_cols, dir_grid, dir_mats)
endif

!****************************************************************
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )
!%===================================
!% solve for forced internal tides
!%===================================
if ( (P%itd_scheme==20)) then
	call forced_itm_solver(trim(P%cpts), P, GD, dir_grid, dir_cols, dir_mats, dir_sols)
!%===================================
else
!%===================================
!% solve for a global barotropic tide
!%===================================
if ((P%fr_scheme <= 2).and.(P%sal_scheme <= 1).and.(P%itd_scheme <= 1)) then
  call baro_solver_linear(trim(P%cpts), P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

elseif ((P%fr_scheme == 3).or.(P%sal_scheme >= 2).or.(P%itd_scheme >= 2)) then
  call baro_solver(trim(P%cpts), P, GD, dir_grid, dir_cols, dir_mats, dir_sols, matlab_dir, save_gridid) ! dir_global

else
      print *, "Inconsistent choice of BBL/SAL/ITD scheme"
      stop
end if
endif
!%===================================
  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

    print *, ""
    call disp ('===>>> Total time spent on solving the system: ' &
			    // trim(ADJUSTL(conv_secs(T2-T1))) // ' CPU, ' &
	            //trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//'  Wall')
    print *, ""
!****************************************************************



!  %===========================================
!  % calculate integrals to check the solutions
!  % do for each component
!  %===========================================
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

	ncpts = len(trim(P%cpts))/2
	allocate(di(ncpts))
	do ccpt = 1, ncpts
		cpt=P%cpts(2*ccpt-1:2*ccpt)
	!	Print the integrals in the very end
		if (P%itd_scheme <= 1) then
			call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, 0, dir_grid, dir_cols, dir_mats, dir_sols, di(ccpt))! parametrized ITD
		else
			call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, 1, dir_grid, dir_cols, dir_mats, dir_sols, di(ccpt)) ! ITD calc with itm
		endif
		call show_domain_integrals(cpt, di(ccpt))
	enddo
	call save_di(trim(P%cpts), di, 'di.txt', dir_sols)

  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

    print *, ""
    call disp ('===>>> Calculation of domain integrals: ' &
			    // trim(ADJUSTL(conv_secs(T2-T1))) // ' CPU, ' &
	            //trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//'  Wall')

!    %============================
!    % clean up temprorary files
!    %============================
!	First argument is clean up level: 0 - none, 1 - remove files in 'temp' dirs, 2 - remove ALL files in dir_mats
	call clean_files(1, dir_cols, dir_grid, dir_mats, dir_sols)

!==========================================================================================

end program
