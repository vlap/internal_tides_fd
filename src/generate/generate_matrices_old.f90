module generate_matrices

     use precisions, only: wp, cwp
     use dispmodule
     use my_trigs
     use my_sparse
     use my_sparse_aggregate
     use save_load
     use generate_grid
	 use control

     contains

!==========================================================================================
subroutine generate_global_matrices(dir_grid, dir_cols, dir_mats, P, GD, H_hg)

     character(len=*) :: dir_grid, dir_cols, dir_mats
     type(params) :: P
     type(grid_dims):: GD

     integer, pointer :: H_hg(:, :)


!%=======================================
!% write H, ta, heq etc. on the columns:
!%=======================================

     call write_global_ccols_col(H_hg, GD%nph, GD%nu, GD%nv, GD%nh, dir_grid, dir_cols)
     deallocate(H_hg) ! don't need it anymore, free space
!%====================================
!% generate and write sparse matrices
!%====================================
     call write_global_cmats_col(GD%nph, GD%nu, GD%nv, GD%nh, P%latP, P%omega, dir_grid, dir_cols, dir_mats)

!     deallocate(ta_ug, ta_vg, ph_ug, ph_vg, H_hg)
!     deallocate(up, vp, hp, iu_ug, iv_vg, ih_hg)

     write(*, '("Global grid generation complete")')

end subroutine
!==========================================================================================

subroutine write_baro_mats(dir_cols, dir_mats, P, nu, nv, nh)

     character(len=*) :: dir_cols, dir_mats
     type(params) :: P
     integer :: nu, nv, nh

     type (triplet) :: h2uddph, h2vddta, u2hddph, v2hddta
     type (triplet) :: u2vf, v2uf
	 real(wp), allocatable:: ta_h(:)
     real(wp), allocatable :: H_u(:), H_v(:)!, H_h(:)
     integer, allocatable :: up(:, :), vp(:, :)

     type (triplet) :: mat
     type (csr) :: mat_csr
     integer, pointer :: bcdiag(:)

!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
! load up, vp, hp
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')

     call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')
     call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
     call load_alloc_vector(H_v, dir_cols // 'H_v.dat')

     call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
     call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
     call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
     call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')
     call load_alloc_sparse(u2vf, dir_mats // 'u2vf.dat')
     call load_alloc_sparse(v2uf, dir_mats // 'v2uf.dat')

!%================================================
!% Collect all the sparse matrices together in mat
!%================================================
      call baro_uvhT_mat(mat, bcdiag, h2uddph, h2vddta, u2hddph, v2hddta, u2vf, v2uf, ta_h, up, vp, nu, nv, nh, nu+nv+nh,&
                         H_u, H_v, P%g, P%re, P%beta0, P%cd, P%Qbar)
!%======================================
!% now store in a temporary matrix file
!%======================================
      call save_sparse(mat, dir_mats // 'temp/mat_init.dat')
      call save_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')

     ! Convert mat from COO to CSR:
     call coo2csr(mat, mat_csr, P%lib)
     call save_sparse(mat_csr, dir_mats // 'temp/mat_csr_init.dat')

     call dealloc_sparse(mat_csr)
     call dealloc_sparse(mat)
     deallocate(bcdiag)

end subroutine

!==========================================================================================
subroutine write_global_cmats_col(nph, nu, nv, nh, latP, omega, dir_grid, dir_cols, dir_mats)
!%===============================================================
!% load the grid: H_hg, ta_ug, ta_vg, ph_ug, ph_vg, th_ug, th_vg,
!%                iu_ug, iv_vg, ih_hg
!%===============================================================
    implicit none

!% Generate some sparse matrices for operations on an Arakawa C grid.
!%
!%   h2uddph, h2u
!%   h2vddta, h2v
!%   u2hddph, u2h
!%   v2hddta, v2h
!%   u2v, u2vf, u2vfsp
!%   v2u, v2uf, v2ufsp

     integer                   ::     nph!, nth
     real(wp), allocatable     ::     th_vg(:), ph_ug(:)

     integer, allocatable :: up(:, :), vp(:, :), hp(:, :)
     integer, intent(in)  ::     nu, nv, nh!, np
     integer              :: nmax
     integer, allocatable :: iu_ug(:, :), iv_vg(:, :), ih_hg(:, :)
     real(wp), allocatable :: H_u(:), H_v(:)!, H_h(:)


     integer, pointer :: ivals(:), jvals(:), ispvals(:), jspvals(:)
     real(wp), pointer :: vals(:), vals1(:), fvals(:), fspvals(:)

     real(wp)  ::      dph, dta, th0
     integer ::      istatus, statusj, status1, status2
     integer ::      cu, cv, ch, cth, cph
     integer ::      cphl, cphr, ctha, cthb, nvs
     integer ::      hpl, hpr, hpa, hpb
     integer ::      upl, upr
     integer ::      vpa, vpb
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     real(wp), intent(in) :: omega, latP

     character(len = *) :: dir_grid, dir_cols, dir_mats

!     real(wp), dimension(4)   :: thc, phc

     type (triplet) :: h2uddph, h2u, h2vddta, h2v, u2hddph, u2h, v2hddta, v2h
     type (triplet) :: u2v, u2vf, u2vfsp, v2u, v2uf, v2ufsp

!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
! load up, vp, hp
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

     call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call load_alloc_vector(th_vg, dir_grid // 'ta_vg.dat') ! actually in Mercator we are using ta but not th
! load iu_ug, iv_vg, ih_hg
     call load_alloc_matrix(iu_ug, dir_grid // 'iu_ug.dat')
     call load_alloc_matrix(iv_vg, dir_grid // 'iv_vg.dat')
     call load_alloc_matrix(ih_hg, dir_grid // 'ih_hg.dat')
! Load columns H_u, H_v
     call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
     call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
!%=================
!% make shortcuts:
!%=================

th0 = d2r(latP)
!ph0 = d2r(lonP)

    dph=2*pi/nph;    ! the distance between ph gridpoints
    dta=dph;    ! the distance between ta gridpoints

!     SECOND ORDER FINITE DIFFERENCE SCHEME

write(*, '("Making matrices:")')

nmax = max(nu, nv, nh)
    if (associated(ivals)) deallocate(ivals) ! 4*nmax for 4th order scheme
    allocate(ivals(2*nmax), stat = istatus)
    if (associated(jvals)) deallocate(jvals)
    allocate(jvals(2*nmax), stat = statusj)
    if (associated(vals)) deallocate(vals)
    allocate(vals(2*nmax), stat = status1)
    if (associated(vals1)) deallocate(vals1)
    allocate(vals1(2*nmax), stat = status2)

    ivals = 0
    jvals = 0
    vals = 0
    vals1 = 0

    nvs = 0


write(*, '(" - d/dph and 1 for h-grid functions, evaluating on the u-grid............. ")', advance = 'no')

!%=============
!% initialize:
!%=============
!!!!!!!!!!!!!!       h2uddph, h2u
!%=============

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cu = 1,nu

       cph=up(cu,1)
       cth=up(cu,2)

      cphr=cph
      cphl=1+modulo((cph-1)-1,nph)
       hpl=ih_hg(cphl,cth)
       hpr=ih_hg(cphr,cth)

       if ((hpl /= 0).and.(hpr /= 0)) then !% we are not at a Western/Eastern boundary

      !% 4th order FD
!              cphll=1+mod((cphl-1)-1,nph)
!              cphrr=1+mod((cphr+1)-1,nph);
!
!              hpll=ih_hg(cth,cphll); hprr=ih_hg(cth,cphrr);
!              if hpll > 0 & hprr > 0 & fdflag == 4
!                nvs=nvs+4;
!                ivals(nvs-3:nvs)=cu;
!                jvals(nvs-3:nvs)=[hpll hpl hpr hprr];
!                vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dph;
!                vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!              else
                nvs=nvs+2
                ivals(nvs-1:nvs)=cu
                jvals(nvs-1:nvs)=(/hpl, hpr/)
                vals(nvs-1:nvs)=(/-1., 1./)/dph
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
       endif

    enddo

call init_sparse(h2uddph,ivals,jvals,vals,nu,nh,nvs)
call init_sparse(h2u,ivals,jvals,vals1,nu,nh,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       h2vddta, h2v
!%=============

write(*, '(" - d/dta and 1 for h grid functions, evaluating on the v-grid............. ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cv = 1,nv

       cph=vp(cv,1)
       cth=vp(cv,2)

      ctha=cth
      cthb=cth-1
       hpa=ih_hg(cph,ctha)
       hpb=ih_hg(cph,cthb)

       if ((hpa /= 0).and.(hpb /= 0)) then !% we are not at a Northern/Southern boundary

      !% 4th order FD
!         cthaa=ctha+1; cthbb=cthb-1;
!         hpaa=ih_hg(cthaa,cph); hpbb=ih_hg(cthbb,cph);
!         if hpaa > 0 && hpbb > 0 && fdflag == 4
!           nvs=nvs+4;
!           ivals(nvs-3:nvs)=cv;
!           jvals(nvs-3:nvs)=[hpbb hpb hpa hpaa];
!           vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dta;
!           vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!         else
                nvs=nvs+2
                ivals(nvs-1:nvs)=cv
                jvals(nvs-1:nvs)=(/hpb, hpa/)
                vals(nvs-1:nvs)=(/-1., 1./)/dta
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
       endif

    enddo

call init_sparse(h2vddta,ivals,jvals,vals,nv,nh,nvs)
call init_sparse(h2v,ivals,jvals,vals1,nv,nh,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')


!%=============
!!!!!!!!!!!!!!       u2hddph, u2h
!%=============

write(*, '(" - d/dph and 1 for u-grid functions, evaluating on the h-grid............. ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do ch = 1,nh

       cph=hp(ch,1)
       cth=hp(ch,2)

      cphl=cph
      cphr=1+modulo((cph+1)-1,nph)
       upl=iu_ug(cphl,cth)
       upr=iu_ug(cphr,cth)


!       if ((hpa /= 0).and.(hpb /= 0)) then !% we are not at a Northern/Southern boundary

      !% 4th order FD
!      cphll=1+mod((cphl-1)-1,nph); cphrr=1+mod((cphr+1)-1,nph);
!      upll=iu_ug(cth,cphll); uprr=iu_ug(cth,cphrr);
!       if upll > 0 & uprr > 0 & fdflag == 4
!         nvs=nvs+4;
!         ivals(nvs-3:nvs)=ch;
!         jvals(nvs-3:nvs)=[upll upl upr uprr];
!         vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dph;
!         vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!       else
                nvs=nvs+2
                ivals(nvs-1:nvs)=ch
                jvals(nvs-1:nvs)=(/upl, upr/)
                vals(nvs-1:nvs)=(/-1., 1./)/dph
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
!       endif

    enddo

call init_sparse(u2hddph,ivals,jvals,vals,nh,nu,nvs)
call init_sparse(u2h,ivals,jvals,vals1,nh,nu,nvs)



call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       v2hddta, v2h
!%=============

write(*, '(" - d/dta and 1 for v-grid functions, evaluating on the h-grid............. ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do ch = 1,nh

       cph=hp(ch,1)
       cth=hp(ch,2)

      ctha=cth+1
      cthb=cth
       vpa=iv_vg(cph,ctha)
       vpb=iv_vg(cph,cthb)

!       cthaa=ctha+1; cthbb=cthb-1;
!       vpaa=iv_vg(cthaa,cph); vpbb=iv_vg(cthbb,cph);
!       if vpaa > 0 & vpbb > 0 & fdflag == 4
!         nvs=nvs+4;
!         ivals(nvs-3:nvs)=ch;
!         jvals(nvs-3:nvs)=[vpbb vpb vpa vpaa];
!         vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dta;
!         vals1(nvs-3:nvs)=[-1 9 9 -1]/18;
!       else
                nvs=nvs+2
                ivals(nvs-1:nvs)=ch
                jvals(nvs-1:nvs)=(/vpb, vpa/)
                vals(nvs-1:nvs)=(/-1., 1./)/dta
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
!       endif

    enddo

call init_sparse(v2hddta,ivals,jvals,vals,nh,nv,nvs)
call init_sparse(v2h,ivals,jvals,vals1,nh,nv,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!**********************************************************************************************************
!**********************************************************************************************************
    deallocate(ivals, jvals, vals, vals1)
    allocate(ivals(4*nmax), jvals(4*nmax), vals(4*nmax), stat = status1)

    if (associated(fvals)) deallocate(fvals)
    allocate(fvals(4*nmax), stat = status2)

    if (associated(ispvals)) deallocate(ispvals) ! nmax for....
    allocate(ispvals(2*nmax), stat = istatus)
    if (associated(jspvals)) deallocate(jspvals)
    allocate(jspvals(2*nmax), stat = statusj)
    if (associated(fspvals)) deallocate(fspvals)
    allocate(fspvals(2*nmax), stat = status1)

    ivals = 0
    jvals = 0
    vals = 0

    ispvals = 0
    jspvals = 0
    fspvals = 0
!**********************************************************************************************************
!**********************************************************************************************************
!%=============
!!!!!!!!!!!!!!       v2u, v2uf, v2ufsp
!%=============

write(*, '(" - v & fv,           evaluating on the u grid............. ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cu = 1,nu

       cph=up(cu,1)
       cth=up(cu,2)

       if (up(cu,3) == 0) then

       ctha=cth+1
       cthb=cth
       cphr=cph
       cphl=1+modulo((cph-1)-1,nph)

                nvs=nvs+4
                ivals(nvs-3:nvs)=cu
                jvals(nvs-3:nvs)=(/iv_vg(cphr,ctha), iv_vg(cphr,cthb), iv_vg(cphl,ctha), iv_vg(cphl,cthb)/)
                vals(nvs-3:nvs)=(/0.25, 0.25, 0.25, 0.25/)

                fvals(nvs-3:nvs)=0.25*calc_f_rot((/th_vg(ctha), th_vg(cthb), th_vg(ctha), th_vg(cthb)/), &
                                                 (/ph_ug(cph), ph_ug(cph), ph_ug(cph), ph_ug(cph)/), &
                                                 th0, omega)

!		    if flags.baro.num.f(2) == 2
!		      cf=1;
!		      while vp(jvals(nvs-4+cf)) == 1
!		        cf=cf+1;
!		      end
!		      ispvals(nvs/4)=cu;
!		      jspvals(nvs/4)=jvals(nvs-4+cf);
!		      fspvals(nvs/4)=4*fvals(nvs-4+cf);
!    elseif flags.baro.num.f(2) == 3
      ispvals(nvs/2-1:nvs/2)=cu;
      jspvals(nvs/2-1:nvs/2)=(/jvals(nvs-3), jvals(nvs)/)
      fspvals(nvs/2-1:nvs/2)=2*(/fvals(nvs-3), fvals(nvs)/)
!    end
       endif

    enddo

call init_sparse(v2u,ivals,jvals,vals,nu,nv,nvs)
call init_sparse(v2uf,ivals,jvals,fvals,nu,nv,nvs)

!if flags.baro.num.f(2) == 2
!  ispvals=ispvals(1:nvs/4);
!  jspvals=jspvals(1:nvs/4);
!  fspvals=fspvals(1:nvs/4);
!  v2ufsp=sparse(ispvals,jspvals,fspvals,nu,nv);
!  save(matfile,'-append','v2ufsp');
!  clear v2ufsp
!elseif flags.baro.num.f(2) == 3

call init_sparse(v2ufsp,ispvals,jspvals,fspvals,nu,nv,nvs/2) ! nvs/2 elements!

!  save(matfile,'-append','v2ufsp');
!  clear v2ufsp
!end


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       u2v, u2vf, u2vfsp
!%=============

write(*, '(" - u & fu,           evaluating on the v-grid............. ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cv = 1,nv

       cph=vp(cv,1)
       cth=vp(cv,2)

       if (vp(cv,3) == 0) then

       ctha=cth
       cthb=cth-1
       cphr=1+modulo((cph+1)-1,nph)
       cphl=cph


                nvs=nvs+4
                ivals(nvs-3:nvs)=cv
                jvals(nvs-3:nvs)=(/iu_ug(cphr,ctha), iu_ug(cphr,cthb), iu_ug(cphl,ctha), iu_ug(cphl,cthb)/)
                vals(nvs-3:nvs)=(/0.25, 0.25, 0.25, 0.25/)
                fvals(nvs-3:nvs)=0.25*calc_f_rot((/th_vg(cth), th_vg(cth), th_vg(cth), th_vg(cth)/), &
                                                 (/ph_ug(cphr), ph_ug(cphr), ph_ug(cphl), ph_ug(cphl)/), &
                                                 th0, omega)

!    if flags.baro.num.f(2) == 2
!      cf=1;
!      while up(jvals(nvs-4+cf)) == 1
!        cf=cf+1;
!      end
!      ispvals(nvs/4)=cv;
!      jspvals(nvs/4)=jvals(nvs-4+cf);
!      fspvals(nvs/4)=4*fvals(nvs-4+cf);
!              elseif flags.baro.num.f(2) == 3
                ispvals(nvs/2-1:nvs/2)=cv;
                jspvals(nvs/2-1:nvs/2)=(/jvals(nvs-3), jvals(nvs)/)
                fspvals(nvs/2-1:nvs/2)=2*(/fvals(nvs-3), fvals(nvs)/)
!              end
       endif

    enddo

call init_sparse(u2v,ivals,jvals,vals,nv,nu,nvs)
call init_sparse(u2vf,ivals,jvals,fvals,nv,nu,nvs)

!if flags.baro.num.f(2) == 2
!  ispvals=ispvals(1:nvs/4);
!  jspvals=jspvals(1:nvs/4);
!  fspvals=fspvals(1:nvs/4);
!  u2vfsp=sparse(ispvals,jspvals,fspvals,nv,nu);
!  save(matfile,'-append','u2vfsp');
!  clear u2vfsp
!elseif flags.baro.num.f(2) == 3

call init_sparse(u2vfsp,ispvals,jspvals,fspvals,nv,nu,nvs/2) ! nvs/2 elements!

!  save(matfile,'-append','u2vfsp');
!  clear u2vfsp
!end


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%================================================
!% Deallocate what's possible
!%================================================
      deallocate(iu_ug, iv_vg, ih_hg, ph_ug, th_vg)

!%===================================================================================================================
!% save sparse matrices h2uddph, h2u, h2vddta, h2v, u2hddph, u2h, v2hddta, v2h, u2v, u2vf, u2vfsp, v2u, v2uf, v2ufsp
!%===================================================================================================================
      call save_sparse(h2uddph, dir_mats // 'h2uddph.dat')
      call save_sparse(h2u, dir_mats // 'h2u.dat')
      call save_sparse(h2vddta, dir_mats // 'h2vddta.dat')
      call save_sparse(h2v, dir_mats // 'h2v.dat')
      call save_sparse(u2hddph, dir_mats // 'u2hddph.dat')
      call save_sparse(u2h, dir_mats // 'u2h.dat')
      call save_sparse(v2hddta, dir_mats // 'v2hddta.dat')
      call save_sparse(v2h, dir_mats // 'v2h.dat')
      call save_sparse(u2v, dir_mats // 'u2v.dat')
      call save_sparse(v2u, dir_mats // 'v2u.dat')

      call save_sparse(u2vfsp, dir_mats // 'u2vfsp.dat')
      call save_sparse(v2ufsp, dir_mats // 'v2ufsp.dat')

!%================================================
!% Adjust the Coriolis matrices to conserve energy
!%================================================
      call write_modified_fmats_col(u2vf, v2uf, H_u, H_v, nu, nv)

      call save_sparse(u2vf, dir_mats // 'u2vf.dat')
      call save_sparse(v2uf, dir_mats // 'v2uf.dat')

end subroutine write_global_cmats_col

!==========================================================================================
!==========================================================================================

subroutine write_modified_fmats_col(u2vf, v2uf, H_u, H_v, nu, nv)
!%====================================================================
!% Different when solving in velocity or volume transport formulations
!%====================================================================
    implicit none

     integer, intent(in)   ::     nu, nv
     real(wp), dimension(:) :: H_u(nu), H_v(nv)

     type (triplet)   :: u2vf, v2uf

!     integer ::      cu, cv, cz

!	Formulation in terms of velocity
!if flags.baro.num.uvform(domflag) == 1
!  v2uf=spdiags(1./sqrt(H_u),0,nu,nu)*v2uf*spdiags(sqrt(H_v),0,nv,nv);
!  u2vf=spdiags(1./sqrt(H_v),0,nv,nv)*u2vf*spdiags(sqrt(H_u),0,nu,nu);
!elseif flags.baro.num.uvform(domflag) == 2
!	Formulation in terms of volume transport

!  v2uf=spdiags(sqrt(H_u),0,nu,nu)*v2uf*spdiags(1./sqrt(H_v),0,nv,nv)
   call left_right_mult_diag(v2uf, real( (H_u)**0.5, wp ), nu, real( 1/((H_v)**0.5), wp ), nv)
!  do cz = 1,v2uf%nz
!       cu = v2uf%indi(cz)
!       cv = v2uf%indj(cz)
!       v2uf%vals(cz) = (H_u(cu))**0.5 * v2uf%vals(cz) / (H_v(cv))**0.5
!  end do

!  u2vf=spdiags(sqrt(H_v),0,nv,nv)*u2vf*spdiags(1./sqrt(H_u),0,nu,nu)
   call left_right_mult_diag(u2vf, real( (H_v)**0.5, wp ), nv, real( 1/((H_u)**0.5), wp ), nu)
!  do cz = 1,u2vf%nz
!       cv = u2vf%indi(cz)
!       cu = u2vf%indj(cz)
!       u2vf%vals(cz) = (H_v(cv))**0.5 * u2vf%vals(cz) / (H_u(cu))**0.5
!  end do

!end

end subroutine write_modified_fmats_col

!==========================================================================================
!==========================================================================================

function calc_f_rot(thc,phc,th0,omega)

implicit none

real(wp), intent(in)           :: thc(:), phc(:), th0, omega
real(wp), dimension(size(phc)) :: calc_f_rot

!% calculates the Coriolis parameters at a point on the
!% computational lat-lon grid. The computational grid is
!% shifted by ph0 in lon, and dropped down th0 from the
!% pole (in the direction of ph0).

calc_f_rot = 2*omega*(-sin(th0)*cos(thc)*cos(phc)+cos(th0)*sin(thc))

end function calc_f_rot

!**********************************************************************************************************
!**********************************************************************************************************

subroutine baro_uvhT_mat(mat, bcdiag, h2uddph, h2vddta, u2hddph, v2hddta, u2vf, v2uf, ta_h, up, vp, nu, nv, nh, np, H_u, H_v, &
                         g, re, beta0, cdg, Q)

    implicit none

!% Produces a linear inversion/propagator matrix for a
!% column uvh using the TRANSPORT VELOCITY formulation.
!% Includes LINEAR friction, scalar SAL;
!% but not frequency dependent d/dt or internal tide drag terms.

!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   ta_h : tau on the h-grid points

!     integer, intent(in)     ::     nph!, nth
!     integer, pointer:: H_hg(:, :)
!     real(wp), pointer :: ta_ug(:), ta_vg(:)!, ph_vg(:)
     integer, allocatable :: up(:, :), vp(:, :)
     integer, intent(in) ::     nu, nv, nh, np
     real(wp), intent(in) :: g, re, beta0, cdg, Q

     real(wp), allocatable ::  ta_h(:)
     real(wp), allocatable :: H_u(:), H_v(:)

     integer, pointer :: p_tmp_int(:)
     real(wp), pointer  :: p_tmp_real(:)
     integer, target  :: iu(nu), iv(nv), ih(nh), ip(np)

     integer ::      j
     integer ::      istatus!, statusj, statusval
     !integer ::      cu, cv, cth, cph, cphl, cphr, ctha, cthb
     integer ::      nmat, filled ! points how many elements out of nmat are already filled
!     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     type (triplet) :: tempmat, mat
     integer, pointer :: bcdiag(:)
     type (triplet) :: h2uddph, h2vddta, u2hddph, v2hddta
     type (triplet) :: u2vf, v2uf

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (j), j=1,nv ) /)
     ih = (/ ( (j), j=1,nh ) /)
     ip = (/ ( (j), j=1,np ) /)

!%==============================
!% first the propagator matrix:
!%==============================
write(*, '("Preparing matrices............. ")', advance = 'no')

!print *, v2uf%nz, u2vf%nz, h2uddph%nz, h2vddta%nz, u2hddph%nz, nu, nv

   nmat = v2uf%nz + u2vf%nz + h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz + nu + nv

     call init_sparse_0(mat, np,np, nmat)
     filled = 0

!if flags.baro.num.f(jj) == 1
!  load(matfile,'v2uf');
  ! 4-element f
  call init_sparse_0(tempmat, v2uf%ni, v2uf%nj, v2uf%nz)
  tempmat = scalar_mult(v2uf, real(-1, wp), real(0, wp))! multiply by -1

!elseif flags.baro.num.f(jj) == 2 | flags.baro.num.f(jj) == 3
! load(matfile,'v2ufsp');
! tempmat=-v2ufsp; ! iterative treatment for f
! clear v2ufsp
!end
![ivals,jvals,svals]=find(tempmat);
     call concatenate_sparse(mat, tempmat, 0, nu, filled)
     call dealloc_sparse(v2uf)
     call dealloc_sparse(tempmat)

write(*, '("1, ")', advance = 'no')

!load(matfile,'h2uddph'); ! Mercator coords
!if flags.coor == 1
!  load(colfile,'H_u');
    !tempmat=(g/re)*(1-beta0)*spdiags(H_u,0,nu,nu)*h2uddph;
    tempmat=h2uddph ! initialize
    !call init_sparse_0(tempmat, v2uf%ni, v2uf%nj, v2uf%nz)
    call left_mult_diag(tempmat, (g/re)*(1-beta0)*H_u, nu)

	call save_sparse(tempmat, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/mats/' &
	// 'tempmat.dat')
!elseif flags.coor == 2
!  load(colfile,'H_u','ta_u');
!  tempmat=(g/re)*(1-beta0)*spdiags(H_u./cos(ta_u),0,nu,nu)*h2uddph; ! lat-lon
!  clear ta_u
!end
![ivals,jvals,svals]=find(tempmat);
     call concatenate_sparse(mat, tempmat, 0, nu+nv, filled)
!     call dealloc_sparse(h2uddph)
     call dealloc_sparse(tempmat)

write(*, '("2, ")', advance = 'no')

!if flags.baro.num.f(jj) == 1
!  load(matfile,'u2vf');
  ! 4-element f
  tempmat = u2vf
!elseif flags.baro.num.f(jj) == 2 | flags.baro.num.f(jj) == 3
!  load(matfile,'u2vfsp');
!  tempmat=u2vfsp;
!  clear u2vfsp
!end
![ivals,jvals,svals]=find(tempmat);
     call concatenate_sparse(mat, tempmat, nu, 0, filled)
     call dealloc_sparse(u2vf)
     call dealloc_sparse(tempmat)

write(*, '("3, ")', advance = 'no')

!load(matfile,'h2vddta')
!load(colfile,'H_v')
    tempmat=h2vddta ! initialize
    call left_mult_diag(tempmat, (g/re)*(1-beta0)*H_v, nv)

     call concatenate_sparse(mat, tempmat, nu, nu+nv, filled)
!     call dealloc_sparse(h2vddta)
     call dealloc_sparse(tempmat)

write(*, '("4, ")', advance = 'no')
!
!load(matfile,'u2hddph')
!load(colfile,'ta_h')
!if flags.coor == 1
! Mercator coords
    tempmat=u2hddph ! initialize
    call left_mult_diag(tempmat, (1/re)*cosh(ta_h)**2, nh)
!elseif flags.coor == 2 ! lat-lon
!  tempmat=(1/re)*sparse(1:nh,1:nh,1./cos(ta_h),nh,nh)*u2hddph;
!end
    call concatenate_sparse(mat, tempmat, nu+nv, 0, filled)
     call dealloc_sparse(u2hddph)
     call dealloc_sparse(tempmat)

write(*, '("5, ")', advance = 'no')

!load(matfile,'v2hddta')
!if flags.coor == 1
! Mercator coords
    tempmat=v2hddta ! initialize
    call left_mult_diag(tempmat, (1/re)*cosh(ta_h)**2, nh)
!elseif flags.coor == 2 ! lat-lon
!  load(colfile,'ta_v')
!  tempmat=(1/re)*sparse(1:nh,1:nh,1./cos(ta_h),nh,nh)*v2hddta*sparse(1:nv,1:nv,cos(ta_v),nv,nv);
!  clear ta_v
!end
     call concatenate_sparse(mat, tempmat, nu+nv, nu, filled)
!     deallocate(ta_h)
     call dealloc_sparse(v2hddta)
     call dealloc_sparse(tempmat)

write(*, '("6, ")', advance = 'no')

!%=============================
!% add linear bottom friction:
!%=============================

!if flags.baro.fric.scheme(jj) == 1
!  cdg=flags.baro.fric.cd(jj);
!  udg=flags.baro.fric.ubar(jj);
!  mat=mat+cdg*udg*sparse(iu,iu,1./H_u,np,np);
!  mat=mat+cdg*udg*sparse(iv,iv,1./H_v,np,np);
!elseif flags.baro.fric.scheme(jj) == 2
! 1/H^2 bc of the transport formulation

  p_tmp_int => iu
  allocate(p_tmp_real(nu), stat = istatus)
  p_tmp_real = cdg*Q/((H_u)**2)

!  print *, nu, size(H_u),size(H_v), tempmat%ni, tempmat%nj, tempmat%nz

  call init_sparse(tempmat,p_tmp_int,p_tmp_int, p_tmp_real, nu,nu, nu)
  call concatenate_sparse(mat, tempmat, 0, 0, filled)
  call dealloc_sparse(tempmat)
  nullify(p_tmp_int, p_tmp_real)

  p_tmp_int => iv
  allocate(p_tmp_real(nv), stat = istatus)
  p_tmp_real = cdg*Q/((H_v)**2)

  call init_sparse(tempmat,p_tmp_int,p_tmp_int, p_tmp_real, nv,nv, nv)
  call concatenate_sparse(mat, tempmat, nu, nu, filled)
  call dealloc_sparse(tempmat)
  nullify(p_tmp_int, p_tmp_real)
!end
     !deallocate(H_u, H_v)

     nullify(p_tmp_int, p_tmp_real)

write(*, '("7, ")', advance = 'no')

!%==================================================
!% define the b.c. diagonal matrix:
!%==================================================

!load(colfile,'up','vp','hp');
  allocate(bcdiag(np), stat = istatus)
  bcdiag(:) = (/ 1-up(:,3), 1-vp(:,3), (1, j=1,nh) /)

write(*, '("8")', advance = 'no')

!%======================================
!% now store in a temporary matrix file
!%======================================
!save(mat0file,'mat','bcmat')

print *, " Done."


end subroutine baro_uvhT_mat

!**********************************************************************************************************
!**********************************************************************************************************

subroutine write_global_ccols_col(H_hg, nph, nu, nv, nh, dir_grid, dir_cols)

    implicit none

!% Generate some useful vectors when using an Arakawa C grid:
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points

     integer, intent(in) ::     nph!, nth
     integer, intent(in) ::     nu, nv, nh!, np
     integer, pointer    ::     H_hg(:, :)
     real(wp), allocatable  :: ta_ug(:), ta_vg(:)!, ph_vg(:)
     integer, allocatable   :: up(:, :), vp(:, :), hp(:, :)

     real(wp), pointer :: ta_u(:), ta_v(:), ta_h(:)
     real(wp), pointer :: H_u(:), H_v(:), H_h(:)

     integer ::      statusu, statusv, statush
     integer ::      cu, cv, ch, cth, cph, cphl, cphr, ctha, cthb
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     character(len=*) :: dir_grid, dir_cols

!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')


write(*, '("Making vectors:")')
call tostring_set(rfmt='F12.1')

write(*, '(" - tau,          on the grid-points................... ")', advance = 'no')

!%=============
!% initialize:
!%=============

!     tau
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(ta_u)) deallocate(ta_u)
    allocate(ta_u(nu), stat = statusu)
    if (associated(ta_v)) deallocate(ta_v)
    allocate(ta_v(nv), stat = statusv)
    if (associated(ta_h)) deallocate(ta_h)
    allocate(ta_h(nh), stat = statush)

ta_u=ta_ug(up(:,2)) ! 2: theta-index
!write(*, '(" 2 ")')
!print *, nv, shape(vp), maxval(vp(:,2)), nth, size(ta_vg)
ta_v=ta_vg(vp(:,2))
!write(*, '(" 3 ")')

ta_h=ta_ug(hp(:,2))
!write(*, '(" 4 ")', advance = 'no')

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on u-grid
write(*, '(" - H,            evaluating on the u-grid............. ")', advance = 'no')

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_u)) deallocate(H_u)
    allocate(H_u(nu), stat = statusu)

do cu = 1, nu

      cph=up(cu,1)
      cth=up(cu,2)

      cphr=cph
      cphl=1+modulo((cph-1)-1,nph)

      if (up(cu,3) == 0) then
        H_u(cu) = (H_hg(cphl,cth)+H_hg(cphr,cth))/2
      else
        H_u(cu) = max( H_hg(cphl,cth), H_hg(cphr,cth) )
      endif

enddo

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on v-grid
write(*, '(" - H,            evaluating on the v-grid............. ")', advance = 'no')

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_v)) deallocate(H_v)
    allocate(H_v(nv), stat = statusv)

do cv = 1, nv

      cph=vp(cv,1)
      cth=vp(cv,2)

      cthb=cth-1
      ctha=cth

      if (vp(cv,3) == 0) then
        H_v(cv) = (H_hg(cph,cthb)+H_hg(cph,ctha))/2
      else
        H_v(cv) = max( H_hg(cph,ctha), H_hg(cph,cthb) )
      endif

enddo

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on h-grid
write(*, '(" - H,            evaluating on the h-grid............. ")', advance = 'no')

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_h)) deallocate(H_h)
    allocate(H_h(nh), stat = statush)

do ch = 1, nh

      cph=hp(ch,1)
      cth=hp(ch,2)

      H_h(ch) = H_hg(cph, cth)

enddo

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============================================
!% save columns H_u, H_v, H_h, ta_u, ta_v, ta_h
!%=============================================
     call save_vector(H_u, dir_cols // 'H_u.dat')
     call save_vector(H_v, dir_cols // 'H_v.dat')
     call save_vector(H_h, dir_cols // 'H_h.dat')
     call save_vector(ta_u, dir_cols // 'ta_u.dat')
     call save_vector(ta_v, dir_cols // 'ta_v.dat')
     call save_vector(ta_h, dir_cols // 'ta_h.dat')

end subroutine write_global_ccols_col

!==========================================================================================

end module generate_matrices



