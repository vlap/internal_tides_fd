! --------------------------------------------------------------------
! MODULE MY_TRIGS contains all definitions of Directories & Files
! --------------------------------------------------------------------

MODULE  my_paths

!==========================================================================================

!    IMPORTANT Directories & Files
     character(len=*), parameter :: model = 'internal_fd'
!     type(varying_string), protected :: MY_TIDES_DATA_PATH

	! theon
     character(len=*), parameter :: data_dir = '/scratch/1/users/amtvl/Tides/data/'!trim(MY_TIDES_DATA_PATH)!
!    ! xigor
!     character(len=*), parameter :: data_dir = '/scratch/a/users/amtvl/Tides/data/'

     character(len=*), parameter :: base_dir = data_dir // 'LAG/' // model// '/'
     character(len=*), parameter :: in_dir = data_dir // 'LAG/' // 'baro_fd/' // 'in/'
     character(len=*), parameter :: out_dir = base_dir // 'out/'
     character(17)               :: save_gridid

     character(len = len(out_dir) + 17 + 1 + len('global/grid/'))     :: dir_cols, dir_grid, dir_mats, dir_sols
     character(len = len(out_dir) + 17)     :: dir_global

     character(len=*), parameter :: matlab_dir = '~/workspace/Matlab/tides/baro_fd/'

     character(len=*), parameter :: etopo_dir = data_dir // 'ETOPO/'
     character(len=*), parameter :: etopo_file = etopo_dir // 'ETOPO1_Ice_g_gmt4.grd'! original ETOPO NetCDF topo file
     character(len=*), parameter :: nocs_dir = data_dir // 'NOCS/'
!     character(len=*), parameter :: nocs_file = nocs_dir // 'nocs_etopo2.nc'! original NOCS NetCDF topo file
!	Processed topo files input/output
     character(len=*), parameter :: topo_dir_out = in_dir // 'topo/etopo/'
     character(len=*), parameter :: topo_dir_in = in_dir // 'topo/nocs/'
     character(len=*), parameter :: topo_file =  'topo_rot_2.0min_pole_15_-40_corrected.dat' !'topo_rot_6.0min_pole_15_-40_interp.dat' !'topo_rot_2.0min_pole_15_-40.dat' ! topo_rot_2.0min_pole_10_-30.dat
                                                                                  ! topo_rot_2.0min_pole_15_-40.dat
     character(len=*), parameter :: N_data_dir = in_dir // 'ocean_N/'
     character(len=*), parameter :: tidal_dir = in_dir // 'tidal/'

     character(len=*), parameter :: P_file = 'control_file.txt'	! contains parameters of the problem and the solver
     character(len=*), parameter :: GD_file = 'grid_file.txt'	! contains grid parameters

CONTAINS

!subroutine init()
!	character(len=255) 	str_tmp
!
!
!	CALL get_environment_variable("MY_TIDES_DATA_PATH", str_tmp)
!
!	MY_TIDES_DATA_PATH = trim(str_tmp)
!
!end subroutine

END MODULE  my_paths
