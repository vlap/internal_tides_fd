module control

     use precisions, only: dp, wp
     use spline
     use save_load

     implicit none

!********************************************************************
! DEFAULT parameters of the problem and the solver.
! The defaults are over-written by the values from the control file
!********************************************************************
type params
!********************************************************************
! General specifications:
!********************************************************************
    integer :: isave = 0      	! true -- save to /yyyy_mm_dd__HH_MM
                                ! false -- save to /0000_00_00__00_00
	integer :: cleanup=1		! 0 - none, 1 - remove files in 'temp' dirs, 2 - remove ALL files in dir_mats
    integer :: graphics = 0		! 1 to plot the final solution; 2 -- plus convergence plots (requires matlab22 command)
	integer :: messages=0		! 0 -- none, 1 -- detailed timings, 2 -- integrals every iteration, 3 -- solver stats (a lot)
    character(len=100) :: cpts='m2'	!  m2, k1, s2, o1, p1, n2, mf, k2, mm, q1
    integer :: ncpts=1
!	flags.graphics=1;
!	flags.coor=1; % 1 for Mercator, 2 for lat-lon
    integer :: grid_link=1	! 0 -- copy; 1 -- make a hard link (fast)
    character(len=17) :: gridid='0' !  0, or a fileid for previous grid...
    integer :: load_sol=0			! 0 -- 0-init guess; 1 -- use the solution found at gridid_sol;
									! 2 -- load and use the tpxo solution
    character(len=17) :: gridid_sol='0'		!
    character(len=17) :: tpxo_dir='0'		!
!********************************************************************
! Global grid:
!********************************************************************
     integer :: nph = 720! 1/4 degree: 1440 ! 1/6 degree: 2160 !1/8 degree: 2880
     integer :: coor = 1! 1 for Mercator, 2 for lat-lon
     integer :: load_etopo = 0 ! true -- load the original ETOPO/NOCS file
								! false -- use topo_file created earlier
	 integer :: etopo_res = 5400 ! number of lon grid points that etopo should be interpolated down to
     integer :: bays=1;  ! 0=keep; 1=remove
     integer :: inlandseas=1;   ! 0=keep; 1=remove
!	flags.global.filter.masklim=0.4;  % set to ocean when this proportion is exceeded...

     integer :: hmin_scheme=1;  ! 1=simple min depth; 2=local critera
     real(dp):: hmin=10;   ! option 1: a simple minimum depth
     integer :: gppw=4;    ! option 2: min number of Grid Points Per Wavelength
!********************************************************************
! Domain integrals:
!********************************************************************
     integer   :: di_log=1			!  save calculated
     integer   :: baro_sol=1		!  make a copy of the calculated baro-only solution
     real(dp)  :: sh_depth = 1000	! deeper than sh_depth meters is DEEP, the rest is SHALLOW ocean
     integer   :: coast_flag = 0	! 0 -- DO NOT perform the coast-open split in DI; 1 -- calculate coastal areas and integrals
     real(dp)  :: coast_dist = 3e5	! closer to the coast than coast_dist is COASTAL, the rest is OPEN ocean
!********************************************************************
!  Bottom friction scheme:
!********************************************************************
     integer   :: fr_scheme=2		! % 1=linear (cd, ubar), 2=linear (cd, Qbar)
									! % 3=quadratic (cd, Q); 4=local area models
     real(dp)  :: cd=0.0025			! nondimensional drag coefficient
     real(dp)  :: ubar=1			! option 1: average value for friction velocity
     real(dp)  :: Qbar=100			! option 2: average value for transport Q
     integer   :: ndays = 0			! length of projection (in days) (optional, if =0 then it is calculated with calc_projection_period)
     integer   :: nppp = 0			! min number of points per period (optional, if =0 then default ppp = 8 is used)
!********************************************************************
!  SAL scheme:
!********************************************************************
     integer	:: sal_scheme=1	! 0=none, 1=scalar approx, >=2 iteration
                          		! 2:real(beta0), 3: betamin < real(beta) < betamax
     real(dp)	:: beta0 = 0.085! initial value for beta
     real(dp)	:: betamin = 0	! minimum value for beta (iteration)
     real(dp)	:: betamax = 0.3! maximum value for beta (iteration)
     integer	:: ntrunc = 360	! % max number of phi gridpoints for spherical harmonics
     integer	:: save_sht = 0	! 0 - keep the SH polynomials in memory between the iterations
								! 1 - save to disk (saves RAM if ntrunc is large)
     real(dp)	:: sal_avg = 10	! % average over this many degrees when calculating beta
!********************************************************************
! BAROTROPIC Internal Tide Ddrag scheme
!********************************************************************
	integer		:: itd_scheme=0		! 0=none; 1=parameterized; 2=local area modelling
	integer		:: N_form = 1		! 1 for uniform, 2 for Ns/(1-z/NL)
    real(dp)	:: Ns = .02		! optional surface stratification
    real(dp)	:: Nl = 500		! optional stratification lengthscale
    real(dp)	:: itd_coeff = 1		! maximum value for the ITD Q parameter
	integer		:: trapped = 1		! 0=all lats; 1=below crit lat
    integer		:: sht_smooth_H=720	! % truncation order for the spherical harmonic representation of H
	integer		:: smooth_type=1	! % truncation order for the spherical harmonic representation of H
	integer		:: baro_on_smoothed=0 !% 0 - baro problem is solved on usmoothed topo (large variations in shallow regions where H_hg varies); 1 - baro problem is solved on same smoothed grid as ITM

    character(len=50) :: N_data='woa05_1deg_pole_15_-40'	! file for statification data
	!N_data=0! 0 for no data...
	integer		:: N_hor_var=0		!% 1-horizontally averaged stratification for IT modelling
									!% 0-horizontally varying as given in N_data
    real(dp)	:: N_hmin = 300		!% The interval for depths on which the eigenvalue problem is solved, i.e.
    real(dp)	:: N_hmax = 11000		!% solved on [0; H_b] where H_b in [N_hmin; N_hmax]
!************** itd_scheme=2, linear global IT scheme *******************
	integer		:: n_modes = 1		! # of IT modes solved for (cannot be large than # of N data points)
	integer		:: eig_order = 4		! discretization order for the eigenvalue problem
	integer		:: n_cheb = 32		! # of Chebyshev nodes used for approximating cn(H)
	real(dp)	:: bulk_damp=0.05	! In domain: IT damping rate (relative to omega~1e-4 Hz). For example it_damp=1/(4*pi)~0.08 means e-folding after travelling 2 wavelength

	integer		:: igw_gppw=4		! Sponge Layer: min number of Grid Points Per Wavelength for IGWs
	integer		:: cn_scheme=1		! Sponge Layer: cn is calculated with: 1 -- WKB; 2 -- actual cn(H)
	real(dp)	:: sponge_damp=1		! Sponge Layer: internal tide damping rate
	integer		:: sponge_forcing=0	! Sponge Layer: flag whether IT should be forced in the sponge
	integer		:: sponge_smooth=2		! Sponge Layer: 1 -- fill in with nearest neighbour at the boundary; 2 -- relaxation smoothing (solve laplace eq within sponge)
	integer		:: ctw_fix=1		! TEMPRORARY, INVESTIGATE! Normal sponge regions for costally trapped waves (K1) may result in regions of negative conversion (hence, bad convergence)
									! Increasing sponge_damp by a factor of 10 removes those. Hence, an ad-hoc fix has been implemented, that multiplies the damping coefficient by 10 for CTW
	integer		:: smooth_forcing=0	! 0 -- no smoothing; 1 -- smooth (u0,v0) using running average over itm_avg degrees (see below)
!************** itd_scheme=3, mixed global IT scheme *******************
	real(dp)	:: eps_rough=0.01	! The criticality cut-off. Regions with |grad(H)|/tan(alpha)>eps_rough are too rough
									! for IT gen there to be adequately resolved with itd_scheme=2. And itd_scheme=1 is applied there
	integer		:: sponge_itd_param=0	! Sponge Layer: 0 -- no parametrised drag in the sponge; 1 -- retain it.
	integer		:: smooth_itd_param=0	! 0 -- no smooth parametrized ITD drag; 1 -- smooth ||grad(H_residual)|| using running average over itm_avg degrees (see below)
!********************************************************************
! LIBRARIES AND THREADING
!********************************************************************
	character   :: lib = 'm' ! 'm' for MKL lib functions; 's' for SparseKit and ad-hoc functions
	integer		:: omp_num_threads=0	! # of processors with shared memory, default = OMP_GET_MAX_THREADS
	integer		:: mpi_num_nodes  =0	! # of nodes for the distributed-memory parallel solver, default = 1
	integer		:: blas_num_threads =0	! 0 - do not use MKL BLAS; else # of threads for MKL BLAS
!********************************************************************
! MATRIX SOLVER
!********************************************************************
	character(len=7) :: solver='pardiso'! umfpack, pardiso, gmres
!******** Iterative algorithms (with preconditioning)    ****************************************
	character(len=7) :: gmres_prec='ilut'		! Select the preconditioning method for GMRES: ilut (not implemented: ilu0)
	integer 	:: gmres_rest_it=3		! the number of the non-restarted FGMRES iterations.
	integer 	:: gmres_tol=6		! specifies residual tolerance: 10^{-gmres_tol} (default is 1.0D-6)
!*************** MKL, PARDISO ****************************************
!**** Direct-Iterative CGS (Conjugate-Gradients Squared) with LU factorization on a first few steps ******
	integer 	:: pardiso_iterative=0	! 0 for direct solver; 1 incurres LU-preconditioned CGS iterations
	integer 	:: pardiso_iter_tol=6	! CGS iterations with astopping tolerance of 10^{-pardiso_iter_tol}
!******** Direct, with parallel LU factorization and Out-of-Core mode ********
	integer 	:: pardiso_symbolic=1	! 0 ('keep none') - incurres symbolic factorization on every iteration
										! 1 ('keep symbolic') - symbolic factorization on the first iteration only
										! (the result is kept in RAM)
	integer 	:: pardiso_ooc=1		! 1 for out-of-core version (RAM require is reduced, slower)
										! 0 for in-core version (default, faster)
	integer 	:: pardiso_max_ram=20*1024	! the total RAM (Mb) that can be used for
										! storing the matrix factors
	integer 	:: pardiso_max_swap=0*1024	! same for swap memory (Mb)
	integer 	:: pardiso_keep_files=0	! 0 -- files saved in OCC will be removed

!********************************************************************
! BAROTROPIC nonlinear solver iterative scheme
!********************************************************************
	integer		:: itn_max=20		! the max alowed number of iterations
	integer		:: cvg_scheme=1			!% 1 for complex dh, 2 for abs(dh)
    real(dp)	:: cvg_dhbar=0.01		!% domain averaged requirement on the residual
    real(dp)	:: cvg_dhmax=0.1		!% pointwise requirement on the residual
    integer		:: p_fric = 1			! 1 for friction; 2 for over-relaxation
    integer		:: p_itd=1 				! 1 for friction; 2 for over-relaxation
    integer		:: itm_method=13		! First digit: 1 -- (Qu, Qv) = (Du/u, Dv/v); 2 -- (Qu, Qv) = (Du/((u*gradH_u + v*gradH_v)*gradH_u), Du/((u*gradH_u + v*gradH_v)*gradH_v));
										! Second digit: 0 -- no smoothing; 1 -- smooth (Qu, Qv) like in SAL, but keep in the rhs the residuals (Du - betau*u, Dv - betav*v);
										! 				3 -- smooth Du using running average over itm_avg degrees
    real(dp)	:: itm_qmax=0.1			! Linear Qu=Du/u can be locally up to 10^20. Keep it within [-itm_qmax; itm_qmax]
    real(dp)	:: itm_avg=0.25			! % average over this many degrees when calculating beta

    real(dp)	:: p_avg = 0.5			! for averaging between iterations; 1 for no averaging (for the barotropic tide)
    real(dp)	:: p_itm_avg = 0.5		! for averaging between iterations; 1 for no averaging (for barotropic tides)
!********************************************************************
! Physical Parameters:
!********************************************************************
     ! In this configuration the North pole is centered on Greenland (30◦W, 80◦N)
     ! and the South pole on Antarctica (30◦W, 80◦S)
     ! The other option of (40◦W, 75◦N) seems to touch water within finer grid
     real(dp) :: latP = 15!10!15 ! Pole of the rotated spherical coordinates
     real(dp) :: lonP = -40!-30!-40

     real(dp)  :: re = 6.371e6                  ! average earth's radius, m
     real(dp)  :: omega = 7.292115e-5                 ! angular velocity, rad/s
     real(dp)  :: g = 9.80665                 ! surface gravity, m/s^2
     real(dp)  :: rhoe = 5515    ! average density of planet, kg/m^2
     real(dp)  :: rhoo = 1030    ! average density of seawater, kg/m^3
!********************************************************************
end type params


     type tide_params
          real(wp) :: amp,lovef,omega0
     end type


    contains
!********************************************************************
subroutine control_file(file_name, P)

  implicit none

  character(len=*) :: file_name
  ! Control file variables
  type(params) :: P

  character(len=100) :: buffer, label
  integer :: pos, pos_
  integer :: fh
  integer :: ios
  integer :: line

fh = 15
ios = 0
line = 0

  open(fh, file=file_name)
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
		!print *, trim(label), trim(buffer)

	        select case (label)
!********************************************************************
! General specifications:
!********************************************************************
	        case ('cpts')
	           read(buffer, *, iostat=ios) P%cpts
!	           print *, 'Read cpts: ', P%cpts
	           P%ncpts = len(trim(P%cpts))/2
	        case ('grid_link')
	           read(buffer, *, iostat=ios) P%grid_link
!	           print *, 'Read grid_link: ', P%grid_link
	        case ('graphics')
	           read(buffer, *, iostat=ios) P%graphics
!	           print *, 'Read graphics: ', P%graphics
	        case ('messages')
	           read(buffer, *, iostat=ios) P%messages
!	           print *, 'Read messages: ', P%messages
	        case ('isave')
	           read(buffer, *, iostat=ios) P%isave
!	           print *, 'Read isave: ', P%isave
	        case ('cleanup')
	           read(buffer, *, iostat=ios) P%cleanup
	           !print *, 'Read cleanup: ', P%cleanup
	        case ('gridid')
	           read(buffer, *, iostat=ios) P%gridid
	           !print *, 'Read gridid: ', P%gridid
	        case ('load_sol')
	           read(buffer, *, iostat=ios) P%load_sol
	           !print *, 'Read load_sol: ', P%load_sol
	        case ('gridid_sol')
	           read(buffer, *, iostat=ios) P%gridid_sol
	           !print *, 'Read gridid_sol: ', P%gridid_sol
	        case ('tpxo_dir')
	           read(buffer, *, iostat=ios) P%tpxo_dir
	           !print *, 'Read tpxo_dir: ', P%tpxo_dir
!********************************************************************
! Global grid:
!********************************************************************
	        case ('nph')
	           read(buffer, *, iostat=ios) P%nph
	           !print *, 'Read nph: ', P%nph
	        case ('coor')
	           read(buffer, *, iostat=ios) P%coor
	           !print *, 'Read coor: ', P%coor
	        case ('load_etopo')
	           read(buffer, *, iostat=ios) P%load_etopo
	           !print *, 'Read load_etopo: ', P%load_etopo
	        case ('etopo_res')
	           read(buffer, *, iostat=ios) P%etopo_res
	           !print *, 'Read etopo_res: ', P%etopo_res
	        case ('bays')
	           read(buffer, *, iostat=ios) P%bays
	           !print *, 'Read bays: ', P%bays
	        case ('inlandseas')
	           read(buffer, *, iostat=ios) P%inlandseas
	           !print *, 'Read inlandseas: ', P%inlandseas
	        case ('hmin_scheme')
	           read(buffer, *, iostat=ios) P%hmin_scheme
	           !print *, 'Read hmin_scheme: ', P%hmin_scheme
	        case ('hmin')
	           read(buffer, *, iostat=ios) P%hmin
	           !print *, 'Read hmin: ', P%hmin
	        case ('gppw')
	           read(buffer, *, iostat=ios) P%gppw
	           !print *, 'Read gppw: ', P%gppw
!********************************************************************
!  Bottom friction scheme:
!********************************************************************
	        case ('fr_scheme')
	           read(buffer, *, iostat=ios) P%fr_scheme
	           !print *, 'Read fr_scheme: ', P%fr_scheme
	        case ('cd')
	           read(buffer, *, iostat=ios) P%cd
	           !print *, 'Read cd: ', P%cd
	        case ('ubar')
	           read(buffer, *, iostat=ios) P%ubar
	           !print *, 'Read ubar: ', P%ubar
	        case ('Qbar')
	           read(buffer, *, iostat=ios) P%Qbar
	           !print *, 'Read Qbar: ', P%Qbar
	        case ('ndays')
	           read(buffer, *, iostat=ios) P%ndays
	           !print *, 'Read ndays: ', P%ndays
	        case ('nppp')
	           read(buffer, *, iostat=ios) P%nppp
	           !print *, 'Read nppp: ', P%nppp
	        case ('di_log')
	           read(buffer, *, iostat=ios) P%di_log
	           !print *, 'Read di_log: ', P%di_log
	        case ('baro_sol')
	           read(buffer, *, iostat=ios) P%baro_sol
	           !print *, 'Read baro_sol: ', P%baro_sol
	        case ('sh_depth')
	           read(buffer, *, iostat=ios) P%sh_depth
	           !print *, 'Read sh_depth: ', P%Qbar
	        case ('coast_flag')
	           read(buffer, *, iostat=ios) P%coast_flag
	           !print *, 'Read coast_flag: ', P%coast_flag
	        case ('coast_dist')
	           read(buffer, *, iostat=ios) P%coast_dist
	           !print *, 'Read coast_dist: ', P%coast_dist
!********************************************************************
!  SAL scheme:
!********************************************************************
	        case ('sal_scheme')
	           read(buffer, *, iostat=ios) P%sal_scheme
	           !print *, 'Read sal_scheme: ', P%sal_scheme
	        case ('beta0')
	           read(buffer, *, iostat=ios) P%beta0
	           !print *, 'Read beta: ', P%beta0
	        case ('betamin')
	           read(buffer, *, iostat=ios) P%betamin
	           !print *, 'Read betamin: ', P%betamin
	        case ('betamax')
	           read(buffer, *, iostat=ios) P%betamax
	           !print *, 'Read betamax: ', P%betamax
	        case ('ntrunc')
	           read(buffer, *, iostat=ios) P%ntrunc
	           !print *, 'Read ntrunc: ', P%ntrunc
	        case ('save_sht')
	           read(buffer, *, iostat=ios) P%save_sht
	           !print *, 'Read save_sht: ', P%save_sht
	        case ('sal_avg')
	           read(buffer, *, iostat=ios) P%sal_avg
	           !print *, 'Read sal_avg: ', P%sal_avg
!********************************************************************
!  ITD scheme:
!********************************************************************
	        case ('itd_scheme')
	           read(buffer, *, iostat=ios) P%itd_scheme
!	           print *, 'Read itd_scheme: ', P%itd_scheme
	        case ('N_data')
	           read(buffer, *, iostat=ios) P%N_data
!	           print *, 'Read N_data: ', P%N_data
	        case ('N_hor_var')
	           read(buffer, *, iostat=ios) P%N_hor_var
!	           print *, 'Read N_hor_var: ', P%N_hor_var
	        case ('N_hmin')
	           read(buffer, *, iostat=ios) P%N_hmin
!	           print *, 'Read N_hmin: ', P%N_hmin
	        case ('N_hmax')
	           read(buffer, *, iostat=ios) P%N_hmax
!	           print *, 'Read N_hmax: ', P%N_hmax

	        case ('N_form')
	           read(buffer, *, iostat=ios) P%N_form
!	           print *, 'Read N_form: ', P%N_form
	        case ('Ns')
	           read(buffer, *, iostat=ios) P%Ns
!	           print *, 'Read N on the surface: ', P%Ns
	        case ('Nl')
	           read(buffer, *, iostat=ios) P%Nl
!	           print *, 'Read stratification lengthscale: ', P%Nl
	        case ('itd_coeff')
	           read(buffer, *, iostat=ios) P%itd_coeff
!	           print *, 'Read a multiplying coeff in the ITD parametrization: ', P%itd_coeff
	        case ('trapped')
	           read(buffer, *, iostat=ios) P%trapped
!	           print *, 'Read trapped at crit lat: ', P%trapped
	        case ('sht_smooth_H')
	           read(buffer, *, iostat=ios) P%sht_smooth_H
!	           print *, 'Read truncation order for SHT of H: ', P%sht_smooth_H
	        case ('smooth_type')
	           read(buffer, *, iostat=ios) P%smooth_type
!	           print *, 'Read smoothing type of H: ', P%smooth_type
	        case ('baro_on_smoothed')
	           read(buffer, *, iostat=ios) P%baro_on_smoothed
!	           print *, 'Read P%baro_on_smoothed: ', P%baro_on_smoothed

	        case ('n_modes')
	           read(buffer, *, iostat=ios) P%n_modes
!	           print *, 'Read # of cn modes: ', P%n_modes
	        case ('eig_order')
	           read(buffer, *, iostat=ios) P%eig_order
!	           print *, 'Read P%eig_order: ', P%eig_order
	        case ('n_cheb')
	           read(buffer, *, iostat=ios) P%n_cheb
!	           print *, 'Read n_cheb: ', P%n_cheb
	        case ('igw_gppw')
	           read(buffer, *, iostat=ios) P%igw_gppw
!	           print *, 'Read igw_gppw: ', P%igw_gppw
	        case ('cn_scheme')
	           read(buffer, *, iostat=ios) P%cn_scheme
!	           print *, 'Read cn_scheme: ', P%cn_scheme
	        case ('bulk_damp')
	           read(buffer, *, iostat=ios) P%bulk_damp
!	           print *, 'Read bulk_damp: ', P%bulk_damp
	        case ('sponge_damp')
	           read(buffer, *, iostat=ios) P%sponge_damp
!	           print *, 'Read sponge_damp: ', P%sponge_damp
	        case ('sponge_forcing')
	           read(buffer, *, iostat=ios) P%sponge_forcing
!	           print *, 'Read sponge_forcing: ', P%sponge_forcing
	        case ('sponge_smooth')
	           read(buffer, *, iostat=ios) P%sponge_smooth
!	           print *, 'Read sponge_smooth: ', P%sponge_smooth
	        case ('ctw_fix')
	           read(buffer, *, iostat=ios) P%ctw_fix
!	           print *, 'Read ctw_fix: ', P%ctw_fix
	        case ('smooth_forcing')
	           read(buffer, *, iostat=ios) P%smooth_forcing
!	           print *, 'Read smooth_forcing: ', P%smooth_forcing
	        case ('sponge_itd_param')
	           read(buffer, *, iostat=ios) P%sponge_itd_param
!	           print *, 'Read sponge_itd_param: ', P%sponge_itd_param
	        case ('smooth_itd_param')
	           read(buffer, *, iostat=ios) P%smooth_itd_param
!	           print *, 'Read smooth_itd_param: ', P%smooth_itd_param
	        case ('eps_rough')
	           read(buffer, *, iostat=ios) P%eps_rough
!	           print *, 'Read eps_rough: ', P%eps_rough
!********************************************************************
!   Nonlinear solver iterative scheme:
!********************************************************************
	        case ('lib')
	           read(buffer, *, iostat=ios) P%lib
!	           print *, 'Read lib: ', P%lib
	        case ('omp_num_threads')
	           read(buffer, *, iostat=ios) P%omp_num_threads
!	           print *, 'Read omp_num_threads: ', P%omp_num_threads
	        case ('mpi_num_nodes')
	           read(buffer, *, iostat=ios) P%mpi_num_nodes
!	           print *, 'Read mpi_num_nodes: ', P%mpi_num_nodes
	        case ('blas_num_threads')
	           read(buffer, *, iostat=ios) P%blas_num_threads
!	           print *, 'Read blas_num_threads: ', P%blas_num_threads
	        case ('solver')
	           read(buffer, *, iostat=ios) P%solver
!	           print *, 'Read solver: ', P%solver
	        case ('gmres_prec')
	           read(buffer, *, iostat=ios) P%gmres_prec
!	           print *, 'Read gmres_prec: ', P%gmres_prec
	        case ('gmres_rest_it')
	           read(buffer, *, iostat=ios) P%gmres_rest_it
!	           print *, 'Read gmres_rest_it: ', P%gmres_rest_it
	        case ('gmres_tol')
	           read(buffer, *, iostat=ios) P%gmres_tol
!	           print *, 'Read gmres_tol: ', P%gmres_tol
	        case ('pardiso_symbolic')
	           read(buffer, *, iostat=ios) P%pardiso_symbolic
!	           print *, 'Read pardiso_symbolic: ', P%pardiso_symbolic
	        case ('pardiso_iterative')
	           read(buffer, *, iostat=ios) P%pardiso_iterative
!	           print *, 'Read pardiso_iterative: ', P%pardiso_iterative
	        case ('pardiso_iter_tol')
	           read(buffer, *, iostat=ios) P%pardiso_iter_tol
!	           print *, 'Read pardiso_iter_tol: ', P%pardiso_iter_tol
	        case ('pardiso_ooc')
	           read(buffer, *, iostat=ios) P%pardiso_ooc
!	           print *, 'Read pardiso_ooc: ', P%pardiso_ooc
	        case ('pardiso_max_ram')
	           read(buffer, *, iostat=ios) P%pardiso_max_ram
!	           print *, 'Read pardiso_max_ram: ', P%pardiso_max_ram
	        case ('pardiso_max_swap')
	           read(buffer, *, iostat=ios) P%pardiso_max_swap
!	           print *, 'Read pardiso_max_swap: ', P%pardiso_max_swap
	        case ('pardiso_keep_files')
	           read(buffer, *, iostat=ios) P%pardiso_keep_files
!	           print *, 'Read pardiso_keep_files: ', P%pardiso_keep_files

	        case ('itn_max')
	           read(buffer, *, iostat=ios) P%itn_max
!	           print *, 'Read itn_max: ', P%itn_max
	        case ('cvg_scheme')
	           read(buffer, *, iostat=ios) P%cvg_scheme
!	           print *, 'Read cvg_scheme: ', P%cvg_scheme
	        case ('p_fric')
	           read(buffer, *, iostat=ios) P%p_fric
!	           print *, 'Read p_fric: ', P%p_fric
	        case ('p_itd')
	           read(buffer, *, iostat=ios) P%p_itd
!	           print *, 'Read p_itd: ', P%p_itd
	        case ('itm_method')
	           read(buffer, *, iostat=ios) P%itm_method
!	           print *, 'Read itm_method: ', P%itm_method
	        case ('itm_qmax')
	           read(buffer, *, iostat=ios) P%itm_qmax
!	           print *, 'Read itm_qmax: ', P%itm_qmax
	        case ('itm_avg')
	           read(buffer, *, iostat=ios) P%itm_avg
!	           print *, 'Read itm_avg: ', P%itm_avg
	        case ('cvg_dhbar')
	           read(buffer, *, iostat=ios) P%cvg_dhbar
!	           print *, 'Read domain averaged requirement on the residual: ', P%cvg_dhbar
	        case ('cvg_dhmax')
	           read(buffer, *, iostat=ios) P%cvg_dhmax
!	           print *, 'Read pointwise requirement on the residual: ', P%cvg_dhmax
	        case ('p_avg')
	           read(buffer, *, iostat=ios) P%p_avg
!	           print *, 'Read averaging between iterations: ', P%p_avg
	        case ('p_itm_avg')
	           read(buffer, *, iostat=ios) P%p_itm_avg
!	           print *, 'Read averaging between iterations(it): ', P%p_itm_avg
!********************************************************************
! Physical Parameters:
!********************************************************************
	        case ('latP')
	           read(buffer, *, iostat=ios) P%latP
	           !print *, 'Read latP: ', P%latP
	        case ('lonP')
	           read(buffer, *, iostat=ios) P%lonP
	           !print *, 'Read lonP: ', P%lonP
	        case ('re')
	           read(buffer, *, iostat=ios) P%re
	           !print *, 'Read re: ', P%re
	        case ('omega')
	           read(buffer, *, iostat=ios) P%omega
	           !print *, 'Read omega: ', P%omega
	        case ('g')
	           read(buffer, *, iostat=ios) P%g
	           !print *, 'Read g: ', P%g
	        case ('rhoe')
	           read(buffer, *, iostat=ios) P%rhoe
	           !print *, 'Read rhoe: ', P%rhoe
	        case ('rhoo')
	           read(buffer, *, iostat=ios) P%rhoo
	           !print *, 'Read rhoo: ', P%rhoo
!********************************************************************
	        case default
	           !print *, 'Skipping invalid label at line', line
	        end select
        end if
     end if
  end do

end subroutine
!**********************************************************************************************************

subroutine control_consistency(P)
!	Quick consistency check

  implicit none
  type(params) :: P

	if ( (P%itd_scheme==3).and.(P%smooth_type==0) ) then
		print *, "*************************************************************************************************"
		print *, "IT scheme with smooth-rough split topo cannot work if smoothing is turned off (smooth_type=0)."
		print *, "Switching it on: smooth_type=3"
		print *, "*************************************************************************************************"
		P%smooth_type=3
	endif

	if ( (P%itd_scheme==20).and.(P%load_sol==0) ) then
		print *, "*************************************************************************************************"
		print *, "Forced internal tides cannot be calculated if no barotropic solution is provided."
		print *, "Attempt to switch on TPXO forcing: load_sol=2"
		print *, "*************************************************************************************************"
		P%load_sol=2
	endif

	if ( (P%itm_method/10==1).and.((P%itm_method-10==0).or.(P%itm_method-10==1).or.(P%itm_method-10==2)) .or. &
		 (P%itm_method/10==2).and.((P%itm_method-20==0).or.(P%itm_method-20==2)) ) then
	else
		print *, "*************************************************************************************************"
		print *, "Choose a valid method of incorporating dynamical ITD into the baro model."
		print *, "Switched to the default (most stable): itm_method=12"
		print *, "*************************************************************************************************"
		P%itm_method=12
	endif

end subroutine

!********************************************************************
! Consistency check when loading previously generated grid
!********************************************************************
logical function cmpr_P_grids(P, P2)

  implicit none

  type(params) :: P, P2
  logical :: cmpr_P_grids

cmpr_P_grids = 	(P%coor == P2%coor ).and.(P%nph == P2%nph ).and. &
	 			(P%load_etopo == P2%load_etopo ).and.&
	 			(P%inlandseas == P2%inlandseas ).and.(P%bays == P2%bays ).and. &
	 			(P%hmin_scheme == P2%hmin_scheme ).and.(P%hmin == P2%hmin ).and. &
	 			(P%lonP == P2%lonP ).and.(P%latP == P2%latP).and.(P%baro_on_smoothed == P2%baro_on_smoothed)

end function

logical function cmpr_P_cn(P, P2)

  implicit none

  type(params) :: P, P2
  logical :: cmpr_P_cn

cmpr_P_cn = 	(P%N_hor_var == P2%N_hor_var ).and.(P%cn_scheme == P2%cn_scheme ).and. &
	 			(P%n_modes <= P2%n_modes ).and.(P%smooth_type == P2%smooth_type ).and.&
	 			(P%eig_order <= P2%eig_order ).and.(P%n_cheb <= P2%n_cheb ).and. &
	 			(P%N_hmin >= P2%N_hmin ).and.(P%N_hmax <= P2%N_hmax )
end function



!**********************************************************************************************************
function get_params(ncpts, cpts)

	implicit none

     integer, intent(in)				:: ncpts
     character(len=2*ncpts), intent(in)   :: cpts

	 integer 			:: ccpt
     character(len=2)	 :: cpt
     type (tide_params) :: get_params(ncpts)

	do ccpt = 1, ncpts
		cpt = cpts(2*ccpt-1:2*ccpt)
		get_params(ccpt) = get_pars(cpt)
	enddo

end function get_params

function get_pars(cpt)

	implicit none

     character(len=2)   :: cpt
     type (tide_params) :: get_pars

!% returns the tidal get_pars%amplitude, love factor, and frequency
!% for the specified tidal component cpt, which is one of
!%
!%  m2, k1, s2, o1, p1, n2, mf, k2, mm, q1

select case(cpt)

       case ('m2')
       get_pars%amp=0.242334
       get_pars%lovef=0.693
       get_pars%omega0=1.405189d-4

       case ('k1')
       get_pars%amp=0.141565
       get_pars%lovef=0.736
       get_pars%omega0=0.7292117d-4

       case ('s2')
       get_pars%amp=0.112743
       get_pars%lovef=0.693
       get_pars%omega0=1.454441d-4

       case ('o1')
       get_pars%amp=0.100661
       get_pars%lovef=0.695
       get_pars%omega0=0.6759774d-4

       case ('p1')
       get_pars%amp=0.046848
       get_pars%lovef=0.706
       get_pars%omega0=0.7252295d-4

       case ('n2')
       get_pars%amp=0.046397
       get_pars%lovef=0.693
       get_pars%omega0=1.378797d-4

       case ('mf')
       get_pars%amp=0.042041
       get_pars%lovef=0.693
       get_pars%omega0=0.053234d-4

       case ('k2')
       get_pars%amp=0.030684
       get_pars%lovef=0.693
       get_pars%omega0=1.458423d-4

       case ('mm')
       get_pars%amp=0.022191
       get_pars%lovef=0.693
       get_pars%omega0=0.026392d-4

       case ('q1')
       get_pars%amp=0.019273
       get_pars%lovef=0.695
       get_pars%omega0=0.6495854d-4

       case ('t2')
       get_pars%amp=0.1
       get_pars%lovef=0.693
       get_pars%omega0=0

       case default
       get_pars%amp=0
       get_pars%lovef=0
       get_pars%omega0=0
       write(*, '("** Need to enter a valid tidal component.")')
       write(*, '("** See function get_pars")')

end select

end function get_pars
!**********************************************************************************************************

subroutine get_load_love(n, nvals, h, k)
! given n, returns the load love numbers h'_n and k'_n.
	implicit none

	integer	:: n, nvals(n)
	real ( kind = 8 ), allocatable	:: nd(:), hd(:), kd(:), ypp(:)
	real ( kind = 8 )				:: ypval, yppval, zero_dp = 0

	real(wp)				:: h(n), k(n)
	real ( kind = 8 )		:: h_tmp(n), k_tmp(n)
	integer	:: d = size( (/ 1, 2, 3, 4, 5, 6, 8, 10, 18, 32, 56, 100, 180, 325, 550, 1000, 1800, 3000, 10000 /) )



	integer :: istat, j


    allocate(nd(d), hd(d), kd(d), ypp(d), stat = istat)
nd = (/ 1, 2, 3, 4, 5, 6, 8, 10, 18, 32, 56, 100, 180, 325, 550, 1000, 1800, 3000, 10000 /)
hd = - (/ 0.290, 1.001, 1.052, 1.053, 1.088, 1.147, 1.291, 1.433, 1.893, 2.379, 2.753, 3.058, 3.474, 4.107, &
		  4.629, 4.906, 4.953, 4.954, 4.956 /)
kd = -(/ 0., 0.615, 0.585, 0.528, 0.516, 0.535, 0.604, 0.682, 0.952, 1.240, 1.402, 1.461, 1.591, 1.928, 2.249,&
		  2.431, 2.465, 2.468, 2.469 /)


!call interp_lagrange ( dim_num, d, real(nd,wp), hd, n, real(nvals,wp), h )
!call interp_lagrange ( dim_num, d, real(nd,wp), kd, n, real(nvals,wp), k )

call spline_cubic_set ( d, nd, hd, 0, zero_dp, 0, zero_dp, ypp ) !0 - the spline should be a quadratic over the first and last intervals
do j = 1,n
	call spline_cubic_val ( d, nd, hd, ypp, real(nvals(j),kind=8), h_tmp(j), ypval, yppval )
end do

call spline_cubic_set ( d, nd, kd, 0,zero_dp, 0,zero_dp, ypp ) !0 - the spline should be a quadratic over the first and last intervals
do j = 1,n
	call spline_cubic_val ( d, nd, kd, ypp, real(nvals(j), kind=8), k_tmp(j), ypval, yppval )
end do
k_tmp = k_tmp/nvals

! back to the variables in our precision
h = h_tmp
k = k_tmp

!call save_vector(nd, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'nd' // '.dat')
!call save_vector(hd, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'hd' // '.dat')
!call save_vector(nvals, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'nvals' // '.dat')
!call save_vector(h, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'hvals' // '.dat')

end subroutine get_load_love


end module
