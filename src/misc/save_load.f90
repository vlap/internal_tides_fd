module save_load

     use precisions, only: sp, dp, wp, cwp, li
     use my_sparse

	 integer, parameter	:: max_chunk_size = (1024**3)  ! the maximum size of a chunk to be written/read in one go (1 gb in bytes)

	 integer, parameter	:: max_i4_size = max_chunk_size/4 ! max length of an integer(4) to be written/read in one go
	 integer, parameter	:: max_i8_size = max_chunk_size/8 ! max length of an integer(8) to be written/read in one go
	 integer, parameter	:: max_r4_size = max_chunk_size/4 ! max length of an real(4) to be written/read in one go
	 integer, parameter	:: max_r8_size = max_chunk_size/8 ! max length of an real(8) to be written/read in one go

!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!---------------------------- WARNING! ---------------------------------------
!	FUNCTIONS SAVE-LOAD WORKING WITH VECTORS USE LONGINTEGER FOR DIMENSIONS
!-----------------------------------------------------------------------------
!---------------------------- SAVE TOPO ---------------------------------------
     interface save_topo
          module procedure save_topo_real_4
          module procedure save_topo_real_8
          module procedure save_topo_int
          module procedure save_topo_cmplx
     end interface
!---------------------------- SAVE SPARSE ---------------------------------------
     interface save_sparse
          module procedure save_cmat
          module procedure save_cmat_dp
          module procedure save_cmat_int
          module procedure save_cmat_csr
          module procedure save_cmat_cmplx
     end interface  
!---------------------------- SAVE MATRIX ---------------------------------------
     interface save_matrix
          module procedure save_mat_4
          module procedure save_mat_8
          module procedure save_mat_int
          module procedure save_mat_cmplx
     end interface

     interface save_matrix_3d
          module procedure save_mat_3d
          module procedure save_mat_3d_int
     end interface
!---------------------------- SAVE VECTOR ---------------------------------------
     interface save_vector
          module procedure save_vec_4
          module procedure save_vec_8
          module procedure save_vec_int
          module procedure save_vec_cmplx_4
          module procedure save_vec_cmplx_8
     end interface  
!---------------------------- LOAD TOPO ---------------------------------------
     interface load_alloc_topo
          module procedure load_alloc_topo_int
          module procedure load_alloc_topo_real_8
          module procedure load_alloc_topo_real_4
          module procedure load_alloc_topo_cmplx_8
     end interface
     interface load_alloc_topo_3d
          module procedure load_alloc_topo_3d_real_8
          module procedure load_alloc_topo_3d_real_4
     end interface
!---------------------------- LOAD MATRIX ---------------------------------------
     interface load_matrix
          module procedure load_mat
          module procedure load_mat_int
     end interface

     interface load_matrix_3d
          module procedure load_mat_3d
     end interface
!---------------------------- LOAD VECTOR ---------------------------------------
     interface load_vector
          module procedure load_vec_4
          module procedure load_vec_8
          module procedure load_vec_int
          module procedure load_vec_cmplx_4
          module procedure load_vec_cmplx_8
     end interface
!---------------------------- LOAD ALLOC VECTOR ---------------------------------
     interface load_alloc_vector
          module procedure load_alloc_vec_4
          module procedure load_alloc_vec_8
          module procedure load_alloc_vec_int
          module procedure load_alloc_vec_cmplx_4
          module procedure load_alloc_vec_cmplx_8
     end interface
!---------------------------- LOAD ALLOC SPARSE ---------------------------------------
     interface load_alloc_sparse
          module procedure load_alloc_cmat
          module procedure load_alloc_cmat_csr
          module procedure load_alloc_cmat_cmplx
     end interface
!---------------------------- LOAD ALLOC MATRIX ---------------------------------
     interface load_alloc_matrix
          module procedure load_alloc_mat_4
          module procedure load_alloc_mat_8
          module procedure load_alloc_mat_int
     end interface
!---------------------------- LOAD ALLOC MATRIX 3D ---------------------------------
     interface load_alloc_matrix_3d
          module procedure load_alloc_mat_3d
     end interface

!---------------------------- WRITE/READ CHUNKS ---------------------------------------
     interface write_chunks
          module procedure write_chunks_4
          module procedure write_chunks_8
          module procedure write_chunks_int_4
          module procedure write_chunks_int_8
          module procedure write_chunks_cmplx_4
          module procedure write_chunks_cmplx_8
     end interface
     interface read_chunks
          module procedure read_chunks_4
          module procedure read_chunks_8
          module procedure read_chunks_int_4
          module procedure read_chunks_int_8
          module procedure read_chunks_cmplx_4
          module procedure read_chunks_cmplx_8
     end interface

     contains

!====================================================================
!----------------------- UPDATE SOLUTION ----------------------------
!====================================================================
subroutine save_solution(nu, nv, nh, itn, load_sol, cpt, uvh, p_avg, dir_sols)

    implicit none

     integer, intent(in) 		:: nu, nv, nh, itn, load_sol
     complex(cwp), intent(in) 	:: uvh(:)
     character(len=2), intent(in)	:: cpt
     character(len = *), intent(in)	:: dir_sols
     real(dp)					:: p_avg ! for averaging between iterations; 1 is no averaging

     complex(cwp), allocatable	:: u(:), v(:), h(:)
     complex(cwp), allocatable	:: u_old(:), v_old(:), h_old(:)
     integer :: iu(nu), iv(nv), ih(nh), istat, j
!-----------------------

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (nu + j), j=1,nv ) /)
     ih = (/ ( (nu + nv + j), j=1,nh ) /)

	    allocate(u(nu), v(nv), h(nh), stat = istat)
!	    print *, 'nu', nu, 'nv', nv, 'nh', nh
!	    print *, size(uvh)
	    u = uvh(iu)
	    v = uvh(iv)
	    h = uvh(ih)
!	    deallocate(uvh)

	    if (.not.((itn == 1).and.(load_sol==0))) then
			call load_alloc_vector(u_old, dir_sols // cpt // '_u' // '_m0' // '.dat')
			call load_alloc_vector(v_old, dir_sols // cpt // '_v' // '_m0' // '.dat')
			call load_alloc_vector(h_old, dir_sols // cpt // '_h' // '_m0' // '.dat')
		        ! rename solution files to old
	        call system('mv ' // dir_sols // cpt // '_u' // '_m0' // '.dat' // ' ' // dir_sols // cpt // '_u' // '_m0_old' // '.dat')
	        call system('mv ' // dir_sols // cpt // '_v' // '_m0' // '.dat' // ' ' // dir_sols // cpt // '_v' // '_m0_old' // '.dat')
	        call system('mv ' // dir_sols // cpt // '_h' // '_m0' // '.dat' // ' ' // dir_sols // cpt // '_h' // '_m0_old' // '.dat')
	        ! average between the new and old solutions
		      u=p_avg*u+(1-p_avg)*u_old
		      v=p_avg*v+(1-p_avg)*v_old
	    	  h=p_avg*h+(1-p_avg)*h_old
	      	deallocate(u_old, v_old, h_old)
	    endif

	     call save_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
	     call save_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
	     call save_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

	!==========================================================================================
	      deallocate(u, v, h)

end subroutine save_solution

!----------------------- COPY BARO SOLUTION ----------------------------
subroutine copy_baro_solution(cpt, dir_sols)

    implicit none
     character(len=2), intent(in)	:: cpt
     character(len = *), intent(in)	:: dir_sols

        ! rename solution files to old
        call system('cp ' // dir_sols // cpt // '_u' // '_m0' // '.dat' // ' ' // dir_sols // cpt // '_u' // '_m0_baro' // '.dat')
        call system('cp ' // dir_sols // cpt // '_v' // '_m0' // '.dat' // ' ' // dir_sols // cpt // '_v' // '_m0_baro' // '.dat')
        call system('cp ' // dir_sols // cpt // '_h' // '_m0' // '.dat' // ' ' // dir_sols // cpt // '_h' // '_m0_baro' // '.dat')

end subroutine copy_baro_solution

!==========================================================================================

subroutine save_itm_solution(nu, nv, nh, cmode, cpt, uvp, itm_avg, dir_sols)

    implicit none

     integer, intent(in) 		:: nu, nv, nh, cmode
     complex(cwp), intent(in) 	:: uvp(:)
     character(len=2), intent(in)	:: cpt
     character(len = *), intent(in)	:: dir_sols
     real(dp)					:: itm_avg ! for averaging between iterations; 1 is no averaging

     complex(cwp), allocatable	:: u(:), v(:), p(:)
     complex(cwp), allocatable	:: u_old(:), v_old(:), p_old(:)
     integer :: iu(nu), iv(nv), ih(nh), istat, j
     logical					:: prev_itn
!-----------------------

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (nu + j), j=1,nv ) /)
     ih = (/ ( (nu + nv + j), j=1,nh ) /)

	    allocate(u(nu), v(nv), p(nh), stat = istat)
	    u = uvp(iu)
	    v = uvp(iv)
	    p = uvp(ih)
!	    deallocate(uvp)
!
	inquire( file=dir_sols// cpt // '_u' // '_m'//tostring(cmode) // '.dat', exist=prev_itn )
!	    if (itn > 1) then
	    if ( (prev_itn).and.(itm_avg/=1) ) then
			call load_alloc_vector(u_old, dir_sols // cpt // '_u' // '_m'//tostring(cmode) // '.dat')
			call load_alloc_vector(v_old, dir_sols // cpt // '_v' // '_m'//tostring(cmode) // '.dat')
			call load_alloc_vector(p_old, dir_sols // cpt // '_p' // '_m'//tostring(cmode) // '.dat')
		        ! rename solution files to old
	        call system('mv ' // dir_sols // cpt // '_u' // '_m'//tostring(cmode) // '.dat' // ' ' // &
	        					 dir_sols // cpt // '_u' // '_m'//tostring(cmode)//'_old' // '.dat')
	        call system('mv ' // dir_sols // cpt // '_v' // '_m'//tostring(cmode) // '.dat' // ' ' // &
	        					 dir_sols // cpt // '_v' // '_m'//tostring(cmode)//'_old' // '.dat')
	        call system('mv ' // dir_sols // cpt // '_p' // '_m'//tostring(cmode) // '.dat' // ' ' // &
	        					 dir_sols // cpt // '_p' // '_m'//tostring(cmode)//'_old' // '.dat')
	        ! average between the new and old solutions
		      u=itm_avg*u+(1-itm_avg)*u_old
		      v=itm_avg*v+(1-itm_avg)*v_old
	    	  p=itm_avg*p+(1-itm_avg)*p_old
	      	deallocate(u_old, v_old, p_old)
	    endif

	     call save_vector(u, dir_sols // cpt // '_u' // '_m'//tostring(cmode) // '.dat')
	     call save_vector(v, dir_sols // cpt // '_v' // '_m'//tostring(cmode) // '.dat')
	     call save_vector(p, dir_sols // cpt // '_p' // '_m'//tostring(cmode) // '.dat')


	      deallocate(u, v, p)

end subroutine save_itm_solution

!====================================================================
!----------------------- SAVE FOR PLOTTING --------------------------
!====================================================================
subroutine save_sol4plot_cmplx(nph, nth, nvar, varp, var, sol4plot_name)

    implicit none

     integer, intent(in) 		:: nph, nth, nvar
     integer, intent(in) 		:: varp(:, :)
     complex(cwp), intent(in)	:: var(:)

    character(len=*) :: sol4plot_name

      OPEN(10,status='unknown', file=sol4plot_name, form='unformatted', action='write', access='stream')

	  write(10) cwp
      write(10) nph, nth, nvar
      write(10) varp(:,1)
      write(10) varp(:,2)
      write(10) real(var, wp)
      write(10) imag(var)
      CLOSE(10)

end subroutine save_sol4plot_cmplx

subroutine save_sol4plot(nph, nth, nvar, varp, var, sol4plot_name)

    implicit none

     integer, intent(in) 		:: nph, nth, nvar
     integer, intent(in) 		:: varp(:, :)
     real(wp), intent(in)	:: var(:)

    character(len=*) :: sol4plot_name

      OPEN(10,status='unknown', file=sol4plot_name, form='unformatted', action='write', access='stream')

	  write(10) wp
      write(10) nph, nth, nvar
      write(10) varp(:,1)
      write(10) varp(:,2)
      write(10) var
      CLOSE(10)

end subroutine save_sol4plot


!==========================================================================================
!----------------------- SAVE/LOAD C-GRID MATRIX IN TRIPLET FORM --------------------------
!==========================================================================================
subroutine save_cmat(cmat, cmat_file)

    implicit none

    type (triplet)   :: cmat
    character(len=*) :: cmat_file

      OPEN(10,status='unknown', file=cmat_file, form='unformatted', action='write', access='stream')

	  write(10) wp
      write(10) cmat%ni, cmat%nj, cmat%nz
!      write(10) cmat%indi
!      write(10) cmat%indj
!      write(10) cmat%vals
      call write_chunks(cmat%indi, 10)
      call write_chunks(cmat%indj, 10)
      call write_chunks(cmat%vals, 10)

      CLOSE(10)

end subroutine save_cmat

subroutine save_cmat_dp(cmat, cmat_file)

    implicit none

    type (triplet_dp)   :: cmat
    character(len=*) :: cmat_file

      OPEN(10,status='unknown', file=cmat_file, form='unformatted', action='write', access='stream')

	  write(10) dp
      write(10) cmat%ni, cmat%nj, cmat%nz
!      write(10) cmat%indi
!      write(10) cmat%indj
!      write(10) cmat%vals
      call write_chunks(cmat%indi, 10)
      call write_chunks(cmat%indj, 10)
      call write_chunks(cmat%vals, 10)

      CLOSE(10)

end subroutine save_cmat_dp
!********************************************************************************************

subroutine save_cmat_int(cmat, cmat_file)

    implicit none

    type (triplet_int)   :: cmat
    character(len=*) :: cmat_file

      OPEN(10,status='unknown', file=cmat_file, form='unformatted', action='write', access='stream')

      write(10) cmat%ni, cmat%nj, cmat%nz
!      write(10) cmat%indi
!      write(10) cmat%indj
!      write(10) cmat%vals
      call write_chunks(cmat%indi, 10)
      call write_chunks(cmat%indj, 10)
      call write_chunks(cmat%vals, 10)

      CLOSE(10)

end subroutine save_cmat_int

!********************************************************************************************
subroutine save_cmat_csr(cmat, cmat_file)

    implicit none

    type (csr)   :: cmat
    character(len=*) :: cmat_file

      OPEN(10,status='unknown', file=cmat_file, form='unformatted', action='write', access='stream')

	  write(10) wp
      write(10) cmat%ni, cmat%nj, cmat%nz
!      write(10) cmat%indi
!      write(10) cmat%indj
!      write(10) cmat%vals
      call write_chunks(cmat%indi, 10)
      call write_chunks(cmat%indj, 10)
      call write_chunks(cmat%vals, 10)

      CLOSE(10)

end subroutine save_cmat_csr
!********************************************************************************************
subroutine save_cmat_cmplx(cmat, cmat_file)

    implicit none

    type (triplet_cmplx)   :: cmat
    character(len=*) :: cmat_file

      OPEN(10,status='unknown', file=cmat_file, form='unformatted', action='write', access='stream')

	  write(10) cwp
      write(10) cmat%ni, cmat%nj, cmat%nz
!      write(10) cmat%indi
!      write(10) cmat%indj
!      write(10) real(cmat%vals)
!      write(10) imag(cmat%vals)
      call write_chunks(cmat%indi, 10)
      call write_chunks(cmat%indj, 10)
      call write_chunks(cmat%vals, 10)

      CLOSE(10)

end subroutine save_cmat_cmplx
!********************************************************************************************
subroutine load_alloc_cmat(cmat, cmat_file)

    implicit none

    integer          :: ni, nj, nz, p, istat
    real(dp),allocatable :: vals_dp(:)
    real(sp),allocatable :: vals_sp(:)

    type (triplet)   :: cmat
    character(len=*) :: cmat_file

      OPEN(unit=10, file=cmat_file,form='unformatted', access='stream', status='old', action='read')

	  read(10) p
      read(10) ni, nj, nz

          if (.not. allocated(cmat%vals)) then
                 call init_sparse_0(cmat, ni,nj, nz)
          elseif ( (cmat%ni/=ni) .or. (cmat%nj/=nj) .or. (cmat%nz/=nz) ) then
                 call dealloc_sparse(cmat)
                 call init_sparse_0(cmat, ni,nj, nz)
          endif

!      read(10) cmat%indi
!      read(10) cmat%indj
      call read_chunks(cmat%indi, 10)
      call read_chunks(cmat%indj, 10)

      if (p == sp) then
      	allocate(vals_sp(nz), stat=istat)
!      	read(10) vals_sp
      	call read_chunks(vals_sp, 10)
      	cmat%vals = vals_sp
      else
      	allocate(vals_dp(nz), stat=istat)
!      	read(10) vals_dp
      	call read_chunks(vals_dp, 10)
      	cmat%vals = vals_dp
      endif

      cmat%ni = ni
      cmat%nj = nj
      cmat%nz = nz

      CLOSE(10)

end subroutine load_alloc_cmat

!********************************************************************************************
subroutine load_alloc_cmat_csr(cmat, cmat_file)

    implicit none

    integer      :: ni, nj, nz, p
    type (csr)   :: cmat
    character(len=*) :: cmat_file

      OPEN(unit=10, file=cmat_file,form='unformatted', access='stream', status='old', action='read')

	  read(10) p
      read(10) ni, nj, nz

          if (.not. allocated(cmat%vals)) then
                 call init_sparse_0(cmat, ni,nj, nz)
          elseif ( (cmat%ni/=ni) .or. (cmat%nj/=nj) .or. (cmat%nz/=nz) ) then
                 call dealloc_sparse(cmat)
                 call init_sparse_0(cmat, ni,nj, nz)
          endif

!      read(10) cmat%indi
!      read(10) cmat%indj
!      read(10) cmat%vals

      call read_chunks(cmat%indi, 10)
      call read_chunks(cmat%indj, 10)
      call read_chunks(cmat%vals, 10)

      cmat%ni = ni
      cmat%nj = nj
      cmat%nz = nz

      CLOSE(10)

end subroutine load_alloc_cmat_csr
!********************************************************************************************
subroutine load_alloc_cmat_cmplx(cmat_cmplx, cmat_file)

    implicit none

    integer          :: ni, nj, nz, p, istat
    type (triplet_cmplx)   :: cmat_cmplx
    character(len=*) :: cmat_file
!    real(wp),allocatable :: re(:), im(:)

      OPEN(unit=10, file=cmat_file,form='unformatted', access='stream', status='old', action='read')

	  read(10) p
      read(10) ni, nj, nz

          if (.not. allocated(cmat_cmplx%vals)) then
                 call init_sparse_0(cmat_cmplx, ni,nj, nz)
          elseif ( (cmat_cmplx%ni/=ni) .or. (cmat_cmplx%nj/=nj) .or. (cmat_cmplx%nz/=nz) ) then
                 call dealloc_sparse(cmat_cmplx)
                 call init_sparse_0(cmat_cmplx, ni,nj, nz)
          endif
!          allocate(re(nz), im(nz), stat=istat)

!      read(10) cmat_cmplx%indi
!      read(10) cmat_cmplx%indj
!      read(10) re
!      read(10) im
      call read_chunks(cmat_cmplx%indi, 10)
      call read_chunks(cmat_cmplx%indj, 10)
      call read_chunks(cmat_cmplx%vals, 10)

      cmat_cmplx%ni = ni
      cmat_cmplx%nj = nj
      cmat_cmplx%nz = nz

!      cmat_cmplx%vals = cmplx(re, im, kind=cwp)

      CLOSE(10)

end subroutine load_alloc_cmat_cmplx

!==========================================================================================
!---------------------------- SAVE/LOAD MATRIX --------------------------
!==========================================================================================
subroutine save_mat_4(mat, mat_file)

    implicit none

    real(sp), dimension(:,:) :: mat
    character(len=*)       :: mat_file

      OPEN(10,status='unknown', file=mat_file, form='unformatted', action='write', access='stream')

	  write(10) sp
      write(10) shape(mat)
      write(10) mat

      CLOSE(10)

end subroutine save_mat_4
subroutine save_mat_8(mat, mat_file)

    implicit none

    real(dp), dimension(:,:) :: mat
    character(len=*)       :: mat_file

      OPEN(10,status='unknown', file=mat_file, form='unformatted', action='write', access='stream')

	  write(10) dp
      write(10) shape(mat)
      write(10) mat

      CLOSE(10)

end subroutine save_mat_8
!********************************************************************************************
subroutine save_mat_int(mat, mat_file)

    implicit none

    integer, dimension(:,:) :: mat
    character(len=*)       :: mat_file

      OPEN(10,status='unknown', file=mat_file, form='unformatted', action='write', access='stream')

      write(10) shape(mat)
      write(10) mat

      CLOSE(10)

end subroutine save_mat_int
!********************************************************************************************
subroutine save_mat_cmplx(mat, mat_file)

    implicit none

    complex(cwp)		:: mat(:,:)
    character(len=*)	:: mat_file

      OPEN(10,status='unknown', file=mat_file, form='unformatted', action='write', access='stream')

	  write(10) cwp
      write(10) shape(mat)
      write(10) real(mat)
      write(10) imag(mat)

      CLOSE(10)

end subroutine save_mat_cmplx
!********************************************************************************************
subroutine save_mat_3d_int(mat, mat_file)

    implicit none

    integer, dimension(:,:,:)	:: mat
    integer						:: c
    character(len=*)       		:: mat_file

      OPEN(10,status='unknown', file=mat_file, form='unformatted', action='write', access='stream')

!	  write(10) wp
      write(10) shape(mat)
      do c = 1,size(mat, 1)
      	write(10) mat(c, :, :)
      enddo

      CLOSE(10)

end subroutine save_mat_3d_int
!********************************************************************************************
subroutine save_mat_3d(mat, mat_file)

    implicit none

    real(wp), dimension(:,:,:)	:: mat
    integer						:: c
    character(len=*)       		:: mat_file

      OPEN(10,status='unknown', file=mat_file, form='unformatted', action='write', access='stream')

	  write(10) wp
      write(10) shape(mat)
      do c = 1,size(mat, 1)
      	write(10) mat(c, :, :)
      enddo

      CLOSE(10)

end subroutine save_mat_3d
!********************************************************************************************

subroutine load_mat(mat, mat_file, m, n)

    implicit none

    integer, optional        :: m, n
    integer					 :: p
    real(wp), dimension(:,:) :: mat
    character(len=*)         :: mat_file

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) m, n
      read(10) mat

      CLOSE(10)

end subroutine load_mat
!********************************************************************************************
subroutine load_mat_3d(mat, mat_file, lin, min, nin)

    implicit none

    integer            :: l, m, n, p
    integer					:: c
    integer, optional        :: lin, min, nin
    real(wp), dimension(:,:,:) :: mat
    character(len=*)         :: mat_file

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) l, m, n

    if ((present(lin)).and.(.not. l==lin)) then
               print *, "Error in subroutine load_mat_3d, input m doesn't match specified in the file: ", l, "/=", lin
    end if
    if ((present(min)).and.(.not. m==min)) then
               print *, "Error in subroutine load_mat_3d, input m doesn't match specified in the file: ", m, "/=", min
    end if
    if ((present(nin)).and.(.not. n==nin)) then
                   print *, "Error in subroutine load_mat_3d, input n don't match specified in the file: ", n, "/=", nin
    end if

      do c = 1,l
      	read(10) mat(c, :, :)
      enddo

      CLOSE(10)

end subroutine load_mat_3d
!********************************************************************************************

subroutine load_alloc_mat_4(mat, mat_file, min, nin)

    implicit none

    integer            :: istatus, m, n, p
    integer, optional        :: min, nin
    real(sp), allocatable     :: mat(:,:)
    character(len=*)         :: mat_file

	real(dp), allocatable     :: mat_dp(:,:)

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) m, n

    if ((present(min)).and.(.not. m==min)) then
               print *, "Error in subroutine load_alloc_mat, input m doesn't match specified in the file: ", m, "/=", min
    end if
    if ((present(nin)).and.(.not. n==nin)) then
                   print *, "Error in subroutine load_alloc_mat, input n don't match specified in the file: ", n, "/=", nin
    end if

          if (.not. allocated(mat)) then
                 allocate(mat(m,n), stat = istatus)
          elseif ( sum(abs(shape(mat) - (/m,n/))) /= 0 ) then
                 deallocate(mat)
                 allocate(mat(m,n), stat = istatus)
          endif

	if (p == sp) then
		allocate(mat(m,n), stat = istatus)
	elseif (p == dp) then
		allocate(mat_dp(m,n), stat = istatus)
		read(10) mat_dp
		mat = mat_dp
	endif

      CLOSE(10)

end subroutine load_alloc_mat_4
!********************************************************************************************
subroutine load_alloc_mat_8(mat, mat_file, min, nin)

    implicit none

    integer            :: istatus, m, n, p
    integer, optional        :: min, nin
    real(dp), allocatable     :: mat(:,:)
    character(len=*)         :: mat_file

	real(sp), allocatable     :: mat_sp(:,:)

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) m, n

    if ((present(min)).and.(.not. m==min)) then
               print *, "Error in subroutine load_alloc_mat, input m doesn't match specified in the file: ", m, "/=", min
    end if
    if ((present(nin)).and.(.not. n==nin)) then
                   print *, "Error in subroutine load_alloc_mat, input n don't match specified in the file: ", n, "/=", nin
    end if

          if (.not. allocated(mat)) then
                 allocate(mat(m,n), stat = istatus)
          elseif ( sum(abs(shape(mat) - (/m,n/))) /= 0 ) then
                 deallocate(mat)
                 allocate(mat(m,n), stat = istatus)
          endif

	if (p == sp) then
		allocate(mat_sp(m,n), stat = istatus)
		read(10) mat_sp
		mat = mat_sp
	elseif (p == dp) then
		read(10) mat
	endif

      CLOSE(10)

end subroutine load_alloc_mat_8
!********************************************************************************************
subroutine load_alloc_mat_3d(mat, mat_file, lin, min, nin)

    implicit none

    integer            :: istatus, l, m, n, p
    integer, optional        :: lin, min, nin
    real(wp), allocatable     :: mat(:,:,:)
    character(len=*)         :: mat_file

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) l, m, n

    if ((present(lin)).and.(.not. l==lin)) then
               print *, "Error in subroutine load_alloc_mat_3d, input m doesn't match specified in the file: ", l, "/=", lin
    end if
    if ((present(min)).and.(.not. m==min)) then
               print *, "Error in subroutine load_alloc_mat_3d, input m doesn't match specified in the file: ", m, "/=", min
    end if
    if ((present(nin)).and.(.not. n==nin)) then
                   print *, "Error in subroutine load_alloc_mat_3d, input n don't match specified in the file: ", n, "/=", nin
    end if

          if (.not. allocated(mat)) then
                 allocate(mat(l,m,n), stat = istatus)
          elseif ( sum(abs(shape(mat) - (/l,m,n/))) /= 0 ) then
                 deallocate(mat)
                 allocate(mat(l,m,n), stat = istatus)
          endif

      read(10) mat


      CLOSE(10)

end subroutine load_alloc_mat_3d
!********************************************************************************************
subroutine load_mat_int(mat, mat_file, min, nin)

    implicit none

    integer        :: m, n
    integer, optional        :: min, nin
    integer, dimension(:,:)  :: mat
    character(len=*)         :: mat_file

!    integer :: k,l

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) m, n

    if(present(min)) m=min
    if(present(nin)) n=nin

          if (sum(abs(shape(mat) - (/m,n/))) /= 0) then
               print *, "Error in subroutine load_mat_int, shape(mat) doesn't match specified in the file: ", shape(mat), "/=", m,n
                 stop
          endif

!do k = 1, m
!do l = 1, n
!print *, k, l
!read(10) mat(k, l)
!
!enddo
!enddo

      read(10) mat

      CLOSE(10)

end subroutine load_mat_int
!********************************************************************************************
subroutine load_alloc_mat_int(mat, mat_file, min, nin)

    implicit none

    integer            :: istatus, m, n
    integer, optional        :: min, nin
    integer, allocatable     :: mat(:,:)
    character(len=*)         :: mat_file

      OPEN(unit=10, file=mat_file,form='unformatted', access='stream', status='old', action='read')

      read(10) m, n

    if ((present(min)).and.(.not. m==min)) then
               print *, "Error in subroutine load_alloc_mat_int, input m doesn't match specified in the file: ", m, "/=", min
    end if
    if ((present(nin)).and.(.not. n==nin)) then
                   print *, "Error in subroutine load_alloc_mat_int, input n don't match specified in the file: ", n, "/=", nin
    end if

          if (.not. allocated(mat)) then
                 allocate(mat(m,n), stat = istatus)
          elseif ( sum(abs(shape(mat) - (/m,n/))) /= 0 ) then
                 deallocate(mat)
                 allocate(mat(m,n), stat = istatus)
          endif

      read(10) mat


      CLOSE(10)

end subroutine load_alloc_mat_int

!==========================================================================================
!---------------------------- SAVE/LOAD VECTOR --------------------------
!==========================================================================================
subroutine save_vec_4(vec, vec_file)

    implicit none

    real(sp), dimension(:) :: vec
    character(len=*)       :: vec_file

      OPEN(10,status='unknown', file=vec_file, form='unformatted', action='write', access='stream')

      write(10) sp
      write(10) int(size(vec), kind=li)
      write(10) vec

      CLOSE(10)

end subroutine save_vec_4
subroutine save_vec_8(vec, vec_file)

    implicit none

    real(dp), dimension(:) :: vec
    character(len=*)       :: vec_file

      OPEN(10,status='unknown', file=vec_file, form='unformatted', action='write', access='stream')

      write(10) dp
      write(10) int(size(vec), kind=li)
      write(10) vec

      CLOSE(10)

end subroutine save_vec_8
!********************************************************************************************
subroutine save_vec_cmplx_4(vec, vec_file)

    implicit none

    complex(sp), dimension(:) :: vec
    character(len=*)       :: vec_file

      OPEN(10,status='unknown', file=vec_file, form='unformatted', action='write', access='stream')

      write(10) cwp
      write(10) int(size(vec), kind=li)
      write(10) real(vec)
      write(10) imag(vec)

      CLOSE(10)

end subroutine save_vec_cmplx_4
!********************************************************************************************
subroutine save_vec_cmplx_8(vec, vec_file)

    implicit none

    complex(dp), dimension(:) :: vec
    character(len=*)       :: vec_file

      OPEN(10,status='unknown', file=vec_file, form='unformatted', action='write', access='stream')

      write(10) cwp
      write(10) int(size(vec), kind=li)
      write(10) real(vec)
      write(10) imag(vec)

      CLOSE(10)

end subroutine save_vec_cmplx_8
!********************************************************************************************
subroutine save_vec_int(vec, vec_file)

    implicit none

    integer, dimension(:) :: vec
    character(len=*)       :: vec_file

      OPEN(10,status='unknown', file=vec_file, form='unformatted', action='write', access='stream')

      write(10) int(size(vec), kind=li)
      write(10) vec

      CLOSE(10)

end subroutine save_vec_int
!********************************************************************************************
subroutine load_vec_4(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)		 :: n
    integer          :: p, istatus
    real(sp)         :: vec(:)
    character(len=*) :: vec_file
    real(dp), allocatable :: vec_dp(:)

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n
              if(present(nin)) nin=n
		if (p == sp) then
			read(10) vec
		else
			allocate(vec_dp(n), stat = istatus)
			read(10) vec_dp
			vec = vec_dp
		endif

      CLOSE(10)

end subroutine load_vec_4

subroutine load_vec_8(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)		 :: n
    integer          :: p, istatus
    real(dp)         :: vec(:)
    character(len=*) :: vec_file
    real(sp), allocatable :: vec_sp(:)

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n
              if(present(nin)) nin=n

		if (p == dp) then
			read(10) vec
		else
			allocate(vec_sp(n), stat = istatus)
			read(10) vec_sp
			vec = vec_sp
		endif

      CLOSE(10)

end subroutine load_vec_8
!********************************************************************************************
subroutine load_vec_int(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: vec(:)
    character(len=*) :: vec_file

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) n
              if(present(nin)) nin=n

      read(10) vec

      CLOSE(10)

end subroutine load_vec_int
!********************************************************************************************
subroutine load_vec_cmplx_4(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: istatus, p
    real(sp), allocatable  :: vec_r(:), vec_i(:)
    complex(sp)     :: vec(:)
    character(len=*) :: vec_file

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n

    if(present(nin)) nin=n

          if (size(vec) /= n) then
               print *, "Error in subroutine load_vec_cmplx, size(vec) doesn't match specified in the file: ", size(vec), "/=", n
               print *, vec_file
                 stop
          endif

    allocate(vec_r(n),vec_i(n), stat = istatus)

      read(10) vec_r
      read(10) vec_i

      CLOSE(10)

      vec = cmplx(vec_r, vec_i, kind = sp)

end subroutine load_vec_cmplx_4
!********************************************************************************************
subroutine load_vec_cmplx_8(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: istatus, p
    real(dp), allocatable  :: vec_r(:), vec_i(:)
    complex(dp)     :: vec(:)
    character(len=*) :: vec_file

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n

    if(present(nin)) nin=n

          if (size(vec) /= n) then
               print *, "Error in subroutine load_vec_cmplx, size(vec) doesn't match specified in the file: ", size(vec), "/=", n
               print *, vec_file
                 stop
          endif

    allocate(vec_r(n),vec_i(n), stat = istatus)

      read(10) vec_r
      read(10) vec_i

      CLOSE(10)

      vec = cmplx(vec_r, vec_i, kind = cwp)

end subroutine load_vec_cmplx_8
!********************************************************************************************
subroutine load_alloc_vec_4(vec, vec_file, nin)

    implicit none

    integer            :: istatus, p
    integer(li)            :: n
    integer, optional  :: nin
    real(sp), allocatable  :: vec(:)
    character(len=*)       :: vec_file
    real(dp), allocatable :: vec_dp(:)

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n

          if (.not. allocated(vec)) then
                 allocate(vec(n), stat = istatus)
          elseif (size(vec) /= n) then
                 deallocate(vec)
                 allocate(vec(n), stat = istatus)
          endif

		if (p == sp) then
			read(10) vec
		else
			allocate(vec_dp(n), stat = istatus)
			read(10) vec_dp
			vec = vec_dp
		endif

      CLOSE(10)

    if(present(nin)) nin=n

end subroutine load_alloc_vec_4

subroutine load_alloc_vec_8(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: istatus, p
    real(dp), allocatable  :: vec(:)
    character(len=*)       :: vec_file
    real(sp), allocatable :: vec_sp(:)

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n

          if (.not. allocated(vec)) then
                 allocate(vec(n), stat = istatus)
          elseif (size(vec) /= n) then
                 deallocate(vec)
                 allocate(vec(n), stat = istatus)
          endif

		if (p == dp) then
			read(10) vec
		else
			allocate(vec_sp(n), stat = istatus)
			read(10) vec_sp
			vec = vec_sp
		endif

      CLOSE(10)

    if(present(nin)) nin=n

end subroutine load_alloc_vec_8
!********************************************************************************************
subroutine load_alloc_vec_int(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: istatus, p
    integer, allocatable   :: vec(:)
    character(len=*)       :: vec_file

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) n

          if (.not. allocated(vec)) then
                 allocate(vec(n), stat = istatus)
          elseif (size(vec) /= n) then
                 deallocate(vec)
                 allocate(vec(n), stat = istatus)
          endif

      read(10) vec

      CLOSE(10)

    if(present(nin)) nin=n

end subroutine load_alloc_vec_int

!********************************************************************************************
subroutine load_alloc_vec_cmplx_4(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: istatus, p
    complex(sp), allocatable   :: vec(:)
    real(sp), allocatable  :: vec_r(:), vec_i(:)
    character(len=*)       :: vec_file

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n

          if (.not. allocated(vec)) then
                 allocate(vec(n), stat = istatus)
          elseif (size(vec) /= n) then
                 deallocate(vec)
                 allocate(vec(n), stat = istatus)
          endif

    allocate(vec_r(n),vec_i(n), stat = istatus)

      read(10) vec_r
      read(10) vec_i

      CLOSE(10)

    vec = cmplx(vec_r, vec_i, kind = sp)

    if(present(nin)) nin=n

end subroutine load_alloc_vec_cmplx_4
!********************************************************************************************
subroutine load_alloc_vec_cmplx_8(vec, vec_file, nin)

    implicit none

    integer, optional:: nin
    integer(li)          :: n
    integer          :: istatus, p
    complex(dp), allocatable   :: vec(:)
    real(dp), allocatable  :: vec_r(:), vec_i(:)
    character(len=*)       :: vec_file

      OPEN(unit=10, file=vec_file,form='unformatted', access='stream', status='old', action='read')

      read(10) p
      read(10) n

          if (.not. allocated(vec)) then
                 allocate(vec(n), stat = istatus)
          elseif (size(vec) /= n) then
                 deallocate(vec)
                 allocate(vec(n), stat = istatus)
          endif

    allocate(vec_r(n),vec_i(n), stat = istatus)

      read(10) vec_r
      read(10) vec_i

      CLOSE(10)

    vec = cmplx(vec_r, vec_i, kind = dp)

    if(present(nin)) nin=n

end subroutine load_alloc_vec_cmplx_8
!==========================================================================================
!---------------------------- SAVE/LOAD INTERPOLATED TOPO --------------------------
!==========================================================================================
subroutine save_topo_int(topo_file, numLons, numLats, xValues, yValues, zValues)

    implicit none

     integer            :: numLons, numLats
     real(dp)  :: xValues(:), yValues(:)
     integer   :: zValues(:, :)
     character(len=*)   :: topo_file

      OPEN(10, status='unknown', file=topo_file,form='unformatted', action='write', access='stream')
      write(10) numLons
      write(10) numLats

      write(10) xValues
      write(10) yValues
      write(10) zValues
      CLOSE(10)

end subroutine save_topo_int

subroutine save_topo_real_4(topo_file, numLons, numLats, xValues, yValues, zValues)

    implicit none

     integer            :: numLons, numLats
     real(sp)  :: xValues(:), yValues(:), zValues(:, :)
     character(len=*)   :: topo_file

      OPEN(10, status='unknown', file=topo_file,form='unformatted', action='write', access='stream')

      write(10) sp
      write(10) numLons
      write(10) numLats

      write(10) xValues
      write(10) yValues
      write(10) zValues
      CLOSE(10)

end subroutine save_topo_real_4

subroutine save_topo_real_8(topo_file, numLons, numLats, xValues, yValues, zValues)

    implicit none

     integer            :: numLons, numLats
     real(dp)  :: xValues(:), yValues(:), zValues(:, :)
     character(len=*)   :: topo_file

      OPEN(10, status='unknown', file=topo_file,form='unformatted', action='write', access='stream')

      write(10) dp
      write(10) numLons
      write(10) numLats

      write(10) xValues
      write(10) yValues
      write(10) zValues
      CLOSE(10)

end subroutine save_topo_real_8

subroutine save_topo_cmplx(topo_file, numLons, numLats, xValues, yValues, zValues)

    implicit none

     integer            :: numLons, numLats
     real(wp)  :: xValues(:), yValues(:)
     complex(cwp)  :: zValues(:, :)
     character(len=*)   :: topo_file

      OPEN(10, status='unknown', file=topo_file,form='unformatted', action='write', access='stream')

      write(10) cwp
      write(10) numLons
      write(10) numLats

      write(10) xValues
      write(10) yValues
      write(10) real(zValues, wp)
      write(10) imag(zValues)
      CLOSE(10)

end subroutine save_topo_cmplx


!********************************************************************************************
subroutine load_alloc_topo_int(topo_file, numLons, numLats, xValues, yValues, zValues, msgs)

    implicit none

	 logical			:: msgs
     integer            :: istatus
     integer            :: numLons, numLats
     real(dp), allocatable  :: xValues(:), yValues(:)
     integer, allocatable   :: zValues(:, :)
     character(len=*)   :: topo_file

if (msgs) then
	write(*, '("Loading topography file: ", a, " ")') topo_file
endif

          OPEN(unit=10, file=topo_file,form='unformatted', action='read', access='stream', status='old')

          read(10) numLons
          read(10) numLats

          if (allocated(ZValues)) deallocate(ZValues)
          allocate(ZValues(numLons, numLats), stat = istatus)
          if (allocated(XValues)) deallocate(XValues)
          allocate(XValues(numLons), stat = istatus)
          if (allocated(YValues)) deallocate(YValues)
          allocate(YValues(numLats), stat = istatus)

          read(10) xValues
          read(10) yValues
          read(10) zValues

          CLOSE(10)

!print *, "done."

end subroutine load_alloc_topo_int

!--------------also is used to read the N_data exported from Matab using function save_topo -----------------!
subroutine load_alloc_topo_real_4(topo_file, numLons, numLats, xValues, yValues, zValues, msgs)

    implicit none

	 logical			:: msgs
     integer            :: istatus, p
     integer            :: numLons, numLats
     real(sp), allocatable  :: xValues(:), yValues(:), zValues(:, :)
     character(len=*)   :: topo_file

     real(dp), allocatable  :: xValues_dp(:), yValues_dp(:), zValues_dp(:, :)

if (msgs) then
	write(*, '("Loading topography file: ", a, " ")') topo_file
endif

          OPEN(unit=10, file=topo_file,form='unformatted', action='read', access='stream', status='old')

          read(10) p
          read(10) numLons
          read(10) numLats

          if (allocated(ZValues)) deallocate(ZValues)
          allocate(ZValues(numLons, numLats), stat = istatus)
          if (allocated(XValues)) deallocate(XValues)
          allocate(XValues(numLons), stat = istatus)
          if (allocated(YValues)) deallocate(YValues)
          allocate(YValues(numLats), stat = istatus)

		if (p == sp) then
          read(10) xValues
          read(10) yValues
          read(10) zValues

		elseif (p == dp) then
          if (allocated(ZValues_dp)) deallocate(ZValues_dp)
          allocate(ZValues_dp(numLons, numLats), stat = istatus)
          if (allocated(XValues_dp)) deallocate(XValues_dp)
          allocate(XValues_dp(numLons), stat = istatus)
          if (allocated(YValues_dp)) deallocate(YValues_dp)
          allocate(YValues_dp(numLats), stat = istatus)
          read(10) xValues_dp
          read(10) yValues_dp
          read(10) zValues_dp

          xValues = xValues_dp
          yValues = yValues_dp
          zValues = zValues_dp
        endif

          CLOSE(10)

!if (msgs) print *, "done."

end subroutine load_alloc_topo_real_4

subroutine load_alloc_topo_real_8(topo_file, numLons, numLats, xValues, yValues, zValues, msgs)

    implicit none

	 logical			:: msgs
     integer            :: istatus, p
     integer            :: numLons, numLats
     real(dp), allocatable  :: xValues(:), yValues(:), zValues(:, :)
     character(len=*)   :: topo_file

     real(sp), allocatable  :: xValues_sp(:), yValues_sp(:), zValues_sp(:, :)

if (msgs) then
	write(*, '("Loading topography file: ", a, " ")') topo_file
endif

          OPEN(unit=10, file=topo_file,form='unformatted', action='read', access='stream', status='old')

          read(10) p
          read(10) numLons
          read(10) numLats

          if (allocated(ZValues)) deallocate(ZValues)
          allocate(ZValues(numLons, numLats), stat = istatus)
          if (allocated(XValues)) deallocate(XValues)
          allocate(XValues(numLons), stat = istatus)
          if (allocated(YValues)) deallocate(YValues)
          allocate(YValues(numLats), stat = istatus)

		if (p == sp) then
          if (allocated(ZValues_sp)) deallocate(ZValues_sp)
          allocate(ZValues_sp(numLons, numLats), stat = istatus)
          if (allocated(XValues_sp)) deallocate(XValues_sp)
          allocate(XValues_sp(numLons), stat = istatus)
          if (allocated(YValues_sp)) deallocate(YValues_sp)
          allocate(YValues_sp(numLats), stat = istatus)
          read(10) xValues_sp
          read(10) yValues_sp
          read(10) zValues_sp

          xValues = xValues_sp
          yValues = yValues_sp
          zValues = zValues_sp
		elseif (p == dp) then
          read(10) xValues
          read(10) yValues
          read(10) zValues

        endif

          CLOSE(10)

!if (msgs) print *, "done."

end subroutine load_alloc_topo_real_8

subroutine load_alloc_topo_cmplx_8(topo_file, numLons, numLats, xValues, yValues, zValues, msgs)


    implicit none

	 logical			:: msgs
     integer            :: istatus, p
     integer            :: numLons, numLats
     real(dp), allocatable  :: xValues(:), yValues(:), zValues_r(:, :), zValues_i(:, :)
     complex(dp), allocatable  :: zValues(:, :)
     character(len=*)   :: topo_file

     real(sp), allocatable  :: xValues_sp(:), yValues_sp(:), zValues_sp_r(:, :), zValues_sp_i(:, :)

if (msgs) then
	write(*, '("Loading topography file: ", a, " ")') topo_file
endif

          OPEN(unit=10, file=topo_file,form='unformatted', action='read', access='stream', status='old')

          read(10) p
          read(10) numLons
          read(10) numLats

          if (allocated(ZValues)) deallocate(ZValues)
          allocate(ZValues(numLons, numLats), stat = istatus)
          if (allocated(XValues)) deallocate(XValues)
          allocate(XValues(numLons), stat = istatus)
          if (allocated(YValues)) deallocate(YValues)
          allocate(YValues(numLats), stat = istatus)

		if (p == sp) then
          if (allocated(ZValues_sp_r)) deallocate(ZValues_sp_r)
          allocate(ZValues_sp_r(numLons, numLats), stat = istatus)
          if (allocated(ZValues_sp_i)) deallocate(ZValues_sp_i)
          allocate(ZValues_sp_i(numLons, numLats), stat = istatus)

          if (allocated(XValues_sp)) deallocate(XValues_sp)
          allocate(XValues_sp(numLons), stat = istatus)
          if (allocated(YValues_sp)) deallocate(YValues_sp)
          allocate(YValues_sp(numLats), stat = istatus)
          read(10) xValues_sp
          read(10) yValues_sp
          read(10) ZValues_sp_r
          read(10) ZValues_sp_i

          xValues = xValues_sp
          yValues = yValues_sp
		  zValues = cmplx(zValues_sp_r, zValues_sp_i, kind = dp)
		elseif (p == dp) then
          if (allocated(ZValues_r)) deallocate(ZValues_r)
          allocate(ZValues_r(numLons, numLats), stat = istatus)
          if (allocated(ZValues_i)) deallocate(ZValues_i)
          allocate(ZValues_i(numLons, numLats), stat = istatus)

          read(10) xValues
          read(10) yValues
          read(10) ZValues_r
          read(10) ZValues_i
		  zValues = cmplx(ZValues_r, ZValues_i, kind = dp)
        endif

          CLOSE(10)

end subroutine load_alloc_topo_cmplx_8

!--------------is used to read the N_data exported from Matab using function save_topo_3d -----------------!
subroutine load_alloc_topo_3d_real_4(topo_file, numLons, numLats, numz, xValues, yValues, zValues, topo, msgs)

    implicit none

	 logical			:: msgs
     integer            :: istatus, p
     integer            :: numLons, numLats, numz
     real(sp), allocatable  :: xValues(:), yValues(:), zValues(:), topo(:, :, :)
     character(len=*)   :: topo_file

     real(dp), allocatable  :: xValues_dp(:), yValues_dp(:), zValues_dp(:), topo_dp(:, :, :)

if (msgs) then
	write(*, '("Loading topography file: ", a, " ")') topo_file
endif

          OPEN(unit=10, file=topo_file,form='unformatted', action='read', access='stream', status='old')

          read(10) p
          read(10) numLons
          read(10) numLats
          read(10) numz

          if (allocated(topo)) deallocate(topo)
          allocate(topo(numLons, numLats, numz), stat = istatus)
          if (allocated(ZValues)) deallocate(ZValues)
          allocate(ZValues(numz), stat = istatus)
          if (allocated(XValues)) deallocate(XValues)
          allocate(XValues(numLons), stat = istatus)
          if (allocated(YValues)) deallocate(YValues)
          allocate(YValues(numLats), stat = istatus)

		if (p == sp) then
          read(10) xValues
          read(10) yValues
          read(10) zValues
          read(10) topo

		elseif (p == dp) then
          if (allocated(topo_dp)) deallocate(topo_dp)
          allocate(topo_dp(numLons, numLats, numz), stat = istatus)
          if (allocated(ZValues_dp)) deallocate(ZValues_dp)
          allocate(ZValues_dp(numz), stat = istatus)
          if (allocated(XValues_dp)) deallocate(XValues_dp)
          allocate(XValues_dp(numLons), stat = istatus)
          if (allocated(YValues_dp)) deallocate(YValues_dp)
          allocate(YValues_dp(numLats), stat = istatus)
          read(10) xValues_dp
          read(10) yValues_dp
          read(10) zValues_dp
          read(10) topo_dp

          xValues = xValues_dp
          yValues = yValues_dp
          zValues = zValues_dp
          topo = topo_dp
        endif

          CLOSE(10)

!if (msgs) print *, "done."

end subroutine load_alloc_topo_3d_real_4

subroutine load_alloc_topo_3d_real_8(topo_file, numLons, numLats, numz, xValues, yValues, zValues, topo, msgs)

    implicit none

	 logical			:: msgs
     integer            :: istatus, p
     integer            :: numLons, numLats, numz
     real(dp), allocatable  :: xValues(:), yValues(:), zValues(:), topo(:, :, :)
     character(len=*)   :: topo_file

     real(sp), allocatable  :: xValues_sp(:), yValues_sp(:), zValues_sp(:), topo_sp(:, :, :)

if (msgs) then
	write(*, '("Loading topography file: ", a, " ")') topo_file
endif

          OPEN(unit=10, file=topo_file,form='unformatted', action='read', access='stream', status='old')

          read(10) p
          read(10) numLons
          read(10) numLats
          read(10) numz

          if (allocated(topo)) deallocate(topo)
          allocate(topo(numLons, numLats, numz), stat = istatus)
          if (allocated(ZValues)) deallocate(ZValues)
          allocate(ZValues(numz), stat = istatus)
          if (allocated(XValues)) deallocate(XValues)
          allocate(XValues(numLons), stat = istatus)
          if (allocated(YValues)) deallocate(YValues)
          allocate(YValues(numLats), stat = istatus)

		if (p == sp) then
          if (allocated(topo_sp)) deallocate(topo_sp)
          allocate(topo_sp(numLons, numLats, numz), stat = istatus)
          if (allocated(ZValues_sp)) deallocate(ZValues_sp)
          allocate(ZValues_sp(numz), stat = istatus)
          if (allocated(XValues_sp)) deallocate(XValues_sp)
          allocate(XValues_sp(numLons), stat = istatus)
          if (allocated(YValues_sp)) deallocate(YValues_sp)
          allocate(YValues_sp(numLats), stat = istatus)
          read(10) xValues_sp
          read(10) yValues_sp
          read(10) zValues_sp
          read(10) topo_sp

          xValues = xValues_sp
          yValues = yValues_sp
          zValues = zValues_sp
          topo = topo_sp
		elseif (p == dp) then
          read(10) xValues
          read(10) yValues
          read(10) zValues
          read(10) topo

        endif

          CLOSE(10)

!if (msgs) print *, "done."

end subroutine load_alloc_topo_3d_real_8

!==========================================================================================
!---------------------------- WRITE/READ IN CHUNKS --------------------------
!==========================================================================================
subroutine write_chunks_4(vec, fh)

    implicit none

    real(sp), dimension(:)	:: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_r4_size-1, n)
	    write(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine write_chunks_4
!********************************************************************************************
subroutine write_chunks_8(vec, fh)

    implicit none

    real(dp), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_r8_size-1, n)
	    write(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine write_chunks_8
!********************************************************************************************
subroutine write_chunks_cmplx_4(vec, fh)

    implicit none

    complex(sp), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_r4_size-1, n)
	    write(fh) real(vec(lpos:rpos))
	    write(fh) imag(vec(lpos:rpos))
	    lpos = rpos + 1
	enddo

end subroutine write_chunks_cmplx_4
!********************************************************************************************
subroutine write_chunks_cmplx_8(vec, fh)

    implicit none

    complex(dp), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_r8_size-1, n)
	    write(fh) real(vec(lpos:rpos))
	    write(fh) imag(vec(lpos:rpos))
	    lpos = rpos + 1
	enddo

end subroutine write_chunks_cmplx_8
!********************************************************************************************
subroutine write_chunks_int_4(vec, fh)

    implicit none

    integer, dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_i4_size-1, n)
	    write(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine write_chunks_int_4
!********************************************************************************************
subroutine write_chunks_int_8(vec, fh)

    implicit none

    integer(kind=li), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_i8_size-1, n)
	    write(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine write_chunks_int_8
!********************************************************************************************

subroutine read_chunks_4(vec, fh)

    implicit none

    real(sp), dimension(:)	:: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_r4_size-1, n)
	    read(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine read_chunks_4
!********************************************************************************************
subroutine read_chunks_8(vec, fh)

    implicit none

    real(dp), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_r8_size-1, n)
	    read(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine read_chunks_8
!********************************************************************************************
subroutine read_chunks_cmplx_4(vec, fh)

    implicit none

    complex(sp), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos
	integer					:: istat
    real(sp),allocatable 	:: re(:), im(:)

	n = size(vec)
	lpos = 1
	rpos = 1

    allocate(re(min(n,max_r4_size)), im(min(n,max_r4_size)), stat=istat)

	do while (lpos <= n)
		rpos = min(lpos+max_r4_size-1, n)
	    read(fh) re(1:rpos-lpos+1)
	    read(fh) im(1:rpos-lpos+1)
	    lpos = rpos + 1
	    vec(lpos:rpos) = cmplx(re(1:rpos-lpos+1), im(1:rpos-lpos+1), kind=sp)
	enddo

end subroutine read_chunks_cmplx_4
!********************************************************************************************
subroutine read_chunks_cmplx_8(vec, fh)

    implicit none

    complex(dp), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos
	integer					:: istat
    real(dp),allocatable 	:: re(:), im(:)

	n = size(vec)
	lpos = 1
	rpos = 1

    allocate(re(min(n,max_r8_size)), im(min(n,max_r8_size)), stat=istat)

	do while (lpos <= n)
		rpos = min(lpos+max_r8_size-1, n)
	    read(fh) re(1:rpos-lpos+1)
	    read(fh) im(1:rpos-lpos+1)
	    lpos = rpos + 1
	    vec(lpos:rpos) = cmplx(re(1:rpos-lpos+1), im(1:rpos-lpos+1), kind=dp)
	enddo

end subroutine read_chunks_cmplx_8
!********************************************************************************************
subroutine read_chunks_int_4(vec, fh)

    implicit none

    integer, dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_i4_size-1, n)
	    read(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine read_chunks_int_4
!********************************************************************************************
subroutine read_chunks_int_8(vec, fh)

    implicit none

    integer(kind=li), dimension(:) :: vec
	integer, intent(in)		:: fh
	integer(kind=li)		:: n, lpos, rpos

	n = size(vec)
	lpos = 1
	rpos = 1

	do while (lpos <= n)
		rpos = min(lpos+max_i8_size-1, n)
	    read(fh) vec(lpos:rpos)
	    lpos = rpos + 1
	enddo

end subroutine read_chunks_int_8
!********************************************************************************************
!==========================================================================================
!---------------------------- SAVE DIRECTORIES --------------------------
!==========================================================================================
subroutine make_save_dirs(isave, now, itm_mod, len_out, output_dir, dir_cols, dir_grid, dir_mats, dir_sols, dir_global)

	implicit none

    integer, intent(in)                :: isave
    integer, intent(in)                :: len_out
    character(len=len_out), intent(in) :: output_dir
    logical, intent(in)			:: itm_mod

    character(17)                :: now
    character(len = len_out + 17 + 1)                            :: dir_global
    character(len = len_out + 17 + 1 + len('global/grid/'))     :: dir_cols, dir_grid, dir_mats, dir_sols

if (isave == 1) then
call time_and_date(output_dir, now)
else
now = '0000_00_00__00_00'
endif

dir_global = output_dir // now // '/'

dir_grid = dir_global // 'global/grid/'
dir_cols = dir_global // 'global/cols/'
dir_mats = dir_global // 'global/mats/'
dir_sols = dir_global // 'global/sols/'

if (isave == 1) then
	write(*, '("Saving data as ", a)') now
else
!print *,'rm -rf ' // dir_global
	call system('rm -rf ' // dir_global)
end if
        call system('mkdir -p ' //  dir_grid)
        call system('mkdir -p ' //  dir_grid // 'temp/')
        call system('mkdir -p ' //  dir_cols)
        call system('mkdir -p ' //  dir_cols // 'temp/')
        call system('mkdir -p ' //  dir_mats)
        call system('mkdir -p ' //  dir_mats // 'temp/')
        call system('mkdir -p ' //  dir_sols)
        call system('mkdir -p ' //  dir_sols // 'temp/')

    if (itm_mod) then
        call system('mkdir -p ' //  dir_grid // 'itm/')
!        call system('mkdir -p ' //  dir_grid // 'itm/' // 'temp/')
        call system('mkdir -p ' //  dir_cols // 'itm/')
!        call system('mkdir -p ' //  dir_cols // 'itm/' // 'temp/')
        call system('mkdir -p ' //  dir_mats // 'itm/')
        call system('mkdir -p ' //  dir_mats // 'itm/' // 'temp/')
        call system('mkdir -p ' //  dir_sols // 'itm/')
        call system('mkdir -p ' //  dir_sols // 'itm/' // 'temp/')
    endif

end subroutine make_save_dirs
!**********************************************************************************************************
!====================================================================
!--------------------- CLEAN UP TEMP FILES --------------------------
!====================================================================
subroutine clean_files(level, dir_cols, dir_grid, dir_mats, dir_sols)

	implicit none

    integer, intent(in)          :: level
    character(len=*), intent(in) :: dir_cols, dir_grid, dir_mats, dir_sols

!	Remove files in 'temp' directories
if (level > 0) then
	write(*, '("Cleaning up temporary files ")')
	call system('rm -rf ' // dir_cols // 'temp/')
	call system('rm -rf ' // dir_grid // 'temp/')
	call system('rm -rf ' // dir_mats // 'temp/')
	call system('rm -rf ' // dir_mats // 'itm/temp/')
	call system('rm -rf ' // dir_sols // 'temp/')
	call system('rm -rf ' // dir_sols // 'itm/temp/')
endif

!	Remove dir_mats and ALL files in it
if (level > 1) then
	call system('rm -rf ' // dir_mats)
endif

end subroutine clean_files
!**********************************************************************************************************
subroutine time_and_date(out_dir, now)

	character(len=*), intent(in) :: out_dir
    character(8)  :: date
    character(10) :: time
    character(17) :: now
    logical		  :: dir_e
    integer		  :: mins, delay=1

    ! using keyword arguments
    call date_and_time(DATE=date,TIME=time)
    now = date(1:4) // '_' // date(5:6) // '_' // date(7:8) // '__' // time(1:2) //'_' // time(3:4)
	inquire( file=out_dir//now//'/.', exist=dir_e )

	do while (dir_e)
		read( time(3:4), '(i2.2)' ) mins
		mins = mins + delay
		write( time(3:4), '(i2.2)' ) mins
		now = date(1:4) // '_' // date(5:6) // '_' // date(7:8) // '__' // time(1:2) //'_' // time(3:4)

		inquire( file=out_dir//now//'/.', exist=dir_e )
		delay = delay + 1
	enddo
end subroutine time_and_date
!**********************************************************************************************************
function conv_secs(seconds)

implicit none

  real, intent(in)		:: seconds
  integer				:: seconds_int

  character (len=100)	:: conv_secs
  character (len=100)	:: str

!    integer							:: istat
    integer							:: d, h, m, s, ms !(day, hours, mins, seconds, milliseconds)
	character(len=*), parameter 	:: dstr='d', hstr=':', mstr=':', sstr='.', msstr=''

	seconds_int = seconds
	d = seconds_int/(24*60*60)
	h = modulo(seconds_int, 24*60*60)/(60*60)
	m = modulo(seconds_int, 60*60)/60
	s = modulo(seconds_int, 60)
	ms = floor((seconds-seconds_int)*1000)
!print *, d, h, m,s
	if (d /= 0) then
		write ( str, '(i2, 1x, a, 1x, i2, a, i2.2, a, i2.2)' ) &
	    				d, dstr, h, hstr, m, mstr, s
    elseif (h /= 0) then
		write ( str, '(i2, a, i2.2, a, i2.2)' ) &
	    				h, hstr, m, mstr, s
    elseif (m /= 0) then
		write ( str, '(i2, a, i2.2, a,i1.1)' ) &
	    		m, mstr, s, sstr, ms/100
    elseif (s /= 0) then
		write ( str, '(i2, a,i3.3,a)' ) &
	    		s, sstr, ms, msstr
    else
		write ( str, '(i1, a,i3.3,a)' ) &
	    		s, sstr, ms, msstr
	endif

	conv_secs = TRIM(str)

end function conv_secs
!**********************************************************************************************************
!**********************************************************************************************************

end module save_load
