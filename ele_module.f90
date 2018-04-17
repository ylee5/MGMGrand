

!---------------------------------------*---------------------------------------!
!                       NASA - Goddard Space Flight Center                      !  
!                   Atmospheric Experiment Branch -  Code 699                   ! 
!                                                                               !
!                     MERCURY'S SPACE ENVIRONEMENT SIMULATOR                    !
!                                   MESS-3D                                     !
!                                                                               !
! Written by :Mehdi Benna                                                       !
!                                                                               !
! Version : 2.03 Two-fluids : solar wind ions + solar wind electrons            !
!                                                                               !
! Last modification date : 07-27-2011                                           !
!                                                                               !
!-------------------------------------------------------------------------------!
!                                  ELE_MODULE.F90                               !
!-------------------------------------------------------------------------------!
! Header 
        module ele_module

        use declarationmodule
        use simulationparameters
        use misc_module
        use ieee_arithmetic

        implicit none

        contains 
        
!-------------------------------------------------------------------------------!   
! Subroutine ele_source_normal                                                  !
!-------------------------------------------------------------------------------!  
        subroutine ele_source_normal(lb)

        implicit none 
        integer, intent(in) :: lb

        return
        end subroutine ele_source_normal

!-------------------------------------------------------------------------------!   
! Subroutine ele_source_surface                                                 !
!-------------------------------------------------------------------------------!  
        subroutine ele_source_surface(lb)

        implicit none 
        integer, intent(in) :: lb
     
        return
        end subroutine ele_source_surface

!-------------------------------------------------------------------------------!   
! Subroutine ele_face_normal                                                    !
!-------------------------------------------------------------------------------! 

        subroutine ele_face_normal(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k

        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: ion_nd_face_xm, ion_nd_face_xp, ion_nd_face_ym, ion_nd_face_yp, ion_nd_face_zm, ion_nd_face_zp

! Interpolate centred values to face values -----------------------------------!
        call  ele_interpol_to_faces_normal(lb, istep, G0_ele_prim(:,:,:,:,lb), ele_face_xm(:,:,:,:,lb), ele_face_xp(:,:,:,:,lb), &
              ele_face_ym(:,:,:,:,lb), ele_face_yp(:,:,:,:,lb), ele_face_zm(:,:,:,:,lb), ele_face_zp(:,:,:,:,lb), order)

        if (hall_mhd_flag==.true.) then
           ion_nd_face_xm = (ion_face_xm(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_xm(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_xm(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_xm(1,:,:,:,lb)-ion_face_xm(6,:,:,:,lb)-ion_face_xm(7,:,:,:,lb)-ion_face_xm(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
           ion_nd_face_xp = (ion_face_xp(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_xp(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_xp(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_xp(1,:,:,:,lb)-ion_face_xp(6,:,:,:,lb)-ion_face_xp(7,:,:,:,lb)-ion_face_xp(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
           ion_nd_face_ym = (ion_face_ym(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_ym(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_ym(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_ym(1,:,:,:,lb)-ion_face_ym(6,:,:,:,lb)-ion_face_ym(7,:,:,:,lb)-ion_face_ym(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
           ion_nd_face_yp = (ion_face_yp(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_yp(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_yp(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_yp(1,:,:,:,lb)-ion_face_yp(6,:,:,:,lb)-ion_face_yp(7,:,:,:,lb)-ion_face_yp(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
           ion_nd_face_zm = (ion_face_zm(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_zm(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_zm(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_zm(1,:,:,:,lb)-ion_face_zm(6,:,:,:,lb)-ion_face_zm(7,:,:,:,lb)-ion_face_zm(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
           ion_nd_face_zp = (ion_face_zp(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_zp(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_zp(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_zp(1,:,:,:,lb)-ion_face_zp(6,:,:,:,lb)-ion_face_zp(7,:,:,:,lb)-ion_face_zp(8,:,:,:,lb))*rho_sc/(M_CO2*amu))

           ele_face_xm(2,:,:,:,lb) = ion_face_xm(2,:,:,:,lb)-((J_face_xm(1,:,:,:,lb))/(ion_nd_face_xm*q*a_sc))
           ele_face_xp(2,:,:,:,lb) = ion_face_xp(2,:,:,:,lb)-((J_face_xp(1,:,:,:,lb))/(ion_nd_face_xp*q*a_sc))
           ele_face_ym(2,:,:,:,lb) = ion_face_ym(2,:,:,:,lb)-((J_face_ym(1,:,:,:,lb))/(ion_nd_face_ym*q*a_sc))
           ele_face_yp(2,:,:,:,lb) = ion_face_yp(2,:,:,:,lb)-((J_face_yp(1,:,:,:,lb))/(ion_nd_face_yp*q*a_sc))
           ele_face_zm(2,:,:,:,lb) = ion_face_zm(2,:,:,:,lb)-((J_face_zm(1,:,:,:,lb))/(ion_nd_face_zm*q*a_sc))
           ele_face_zp(2,:,:,:,lb) = ion_face_zp(2,:,:,:,lb)-((J_face_zp(1,:,:,:,lb))/(ion_nd_face_zp*q*a_sc))
           ele_face_xm(3,:,:,:,lb) = ion_face_xm(3,:,:,:,lb)-((J_face_xm(2,:,:,:,lb))/(ion_nd_face_xm*q*a_sc))
           ele_face_xp(3,:,:,:,lb) = ion_face_xp(3,:,:,:,lb)-((J_face_xp(2,:,:,:,lb))/(ion_nd_face_xp*q*a_sc))
           ele_face_ym(3,:,:,:,lb) = ion_face_ym(3,:,:,:,lb)-((J_face_ym(2,:,:,:,lb))/(ion_nd_face_ym*q*a_sc))
           ele_face_yp(3,:,:,:,lb) = ion_face_yp(3,:,:,:,lb)-((J_face_yp(2,:,:,:,lb))/(ion_nd_face_yp*q*a_sc))
           ele_face_zm(3,:,:,:,lb) = ion_face_zm(3,:,:,:,lb)-((J_face_zm(2,:,:,:,lb))/(ion_nd_face_zm*q*a_sc))
           ele_face_zp(3,:,:,:,lb) = ion_face_zp(3,:,:,:,lb)-((J_face_zp(2,:,:,:,lb))/(ion_nd_face_zp*q*a_sc))
           ele_face_xm(4,:,:,:,lb) = ion_face_xm(4,:,:,:,lb)-((J_face_xm(3,:,:,:,lb))/(ion_nd_face_xm*q*a_sc))
           ele_face_xp(4,:,:,:,lb) = ion_face_xp(4,:,:,:,lb)-((J_face_xp(3,:,:,:,lb))/(ion_nd_face_xp*q*a_sc))
           ele_face_ym(4,:,:,:,lb) = ion_face_ym(4,:,:,:,lb)-((J_face_ym(3,:,:,:,lb))/(ion_nd_face_ym*q*a_sc))
           ele_face_yp(4,:,:,:,lb) = ion_face_yp(4,:,:,:,lb)-((J_face_yp(3,:,:,:,lb))/(ion_nd_face_yp*q*a_sc))
           ele_face_zm(4,:,:,:,lb) = ion_face_zm(4,:,:,:,lb)-((J_face_zm(3,:,:,:,lb))/(ion_nd_face_zm*q*a_sc))
           ele_face_zp(4,:,:,:,lb) = ion_face_zp(4,:,:,:,lb)-((J_face_zp(3,:,:,:,lb))/(ion_nd_face_zp*q*a_sc))

        else

           ele_face_xm(2,:,:,:,lb) = ion_face_xm(2,:,:,:,lb)
           ele_face_xp(2,:,:,:,lb) = ion_face_xp(2,:,:,:,lb)
           ele_face_ym(2,:,:,:,lb) = ion_face_ym(2,:,:,:,lb)
           ele_face_yp(2,:,:,:,lb) = ion_face_yp(2,:,:,:,lb)
           ele_face_zm(2,:,:,:,lb) = ion_face_zm(2,:,:,:,lb)
           ele_face_zp(2,:,:,:,lb) = ion_face_zp(2,:,:,:,lb)
           ele_face_xm(3,:,:,:,lb) = ion_face_xm(3,:,:,:,lb)
           ele_face_xp(3,:,:,:,lb) = ion_face_xp(3,:,:,:,lb)
           ele_face_ym(3,:,:,:,lb) = ion_face_ym(3,:,:,:,lb)
           ele_face_yp(3,:,:,:,lb) = ion_face_yp(3,:,:,:,lb)
           ele_face_zm(3,:,:,:,lb) = ion_face_zm(3,:,:,:,lb)
           ele_face_zp(3,:,:,:,lb) = ion_face_zp(3,:,:,:,lb)
           ele_face_xm(4,:,:,:,lb) = ion_face_xm(4,:,:,:,lb)
           ele_face_xp(4,:,:,:,lb) = ion_face_xp(4,:,:,:,lb)
           ele_face_ym(4,:,:,:,lb) = ion_face_ym(4,:,:,:,lb)
           ele_face_yp(4,:,:,:,lb) = ion_face_yp(4,:,:,:,lb)
           ele_face_zm(4,:,:,:,lb) = ion_face_zm(4,:,:,:,lb)
           ele_face_zp(4,:,:,:,lb) = ion_face_zp(4,:,:,:,lb)

        end if

! Convert primitive variables to conservative variables -----------------------!                                                           
        call  ele_primitive_to_conserv(lb, ele_face_xm(:,:,:,:,lb), ele_cons_face_xm(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_xp(:,:,:,:,lb), ele_cons_face_xp(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_ym(:,:,:,:,lb), ele_cons_face_ym(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_yp(:,:,:,:,lb), ele_cons_face_yp(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_zm(:,:,:,:,lb), ele_cons_face_zm(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_zp(:,:,:,:,lb), ele_cons_face_zp(:,:,:,:,lb))

! Compute the average face values ---------------------------------------------!
        ele_face_xm_av(:,:,:,:,lb)=2.
        ele_face_xp_av(:,:,:,:,lb)=2.
        ele_face_ym_av(:,:,:,:,lb)=2.
        ele_face_yp_av(:,:,:,:,lb)=2.
        ele_face_zm_av(:,:,:,:,lb)=2.
        ele_face_zp_av(:,:,:,:,lb)=2.
        ele_face_ex_av(:,:,:,:,lb)=2.

        do i=nguard+1,nxb+nguard
           do j=nguard+1,nyb+nguard
              do k=nguard+1,nzb+nguard
                 ele_face_xm_av(:,i,j,k,lb)=0.5*(ele_face_xm(:,i,j,k,lb)+ele_face_xp(:,i-1,j,k,lb))
                 ele_face_xp_av(:,i,j,k,lb)=0.5*(ele_face_xp(:,i,j,k,lb)+ele_face_xm(:,i+1,j,k,lb))
                 ele_face_ym_av(:,i,j,k,lb)=0.5*(ele_face_ym(:,i,j,k,lb)+ele_face_yp(:,i,j-1,k,lb))
                 ele_face_yp_av(:,i,j,k,lb)=0.5*(ele_face_yp(:,i,j,k,lb)+ele_face_ym(:,i,j+1,k,lb))
                 ele_face_zm_av(:,i,j,k,lb)=0.5*(ele_face_zm(:,i,j,k,lb)+ele_face_zp(:,i,j,k-1,lb))
                 ele_face_zp_av(:,i,j,k,lb)=0.5*(ele_face_zp(:,i,j,k,lb)+ele_face_zm(:,i,j,k+1,lb))
              end do
           end do
        end do

        end subroutine ele_face_normal

!-------------------------------------------------------------------------------!   
! Subroutine ele_face_surface                                                   !
!-------------------------------------------------------------------------------! 

        subroutine ele_face_surface(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k

        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: ion_nd_face_xm, ion_nd_face_xp, ion_nd_face_ym, ion_nd_face_yp, ion_nd_face_zm, ion_nd_face_zp


! Interpolate centred values to face values -----------------------------------!
        call  ele_interpol_to_faces_surface(lb, istep, G0_ele_prim(:,:,:,:,lb), ele_face_xm(:,:,:,:,lb), ele_face_xp(:,:,:,:,lb), &
              ele_face_ym(:,:,:,:,lb), ele_face_yp(:,:,:,:,lb), ele_face_zm(:,:,:,:,lb), ele_face_zp(:,:,:,:,lb), &
              ele_face_ei(:,:,:,:,lb), ele_face_eo(:,:,:,:,lb), order)

        if (hall_mhd_flag==.true.) then

          ion_nd_face_xm = (ion_face_xm(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_xm(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_xm(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_xm(1,:,:,:,lb)-ion_face_xm(6,:,:,:,lb)-ion_face_xm(7,:,:,:,lb)-ion_face_xm(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_xp = (ion_face_xp(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_xp(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_xp(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_xp(1,:,:,:,lb)-ion_face_xp(6,:,:,:,lb)-ion_face_xp(7,:,:,:,lb)-ion_face_xp(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_ym = (ion_face_ym(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_ym(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_ym(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_ym(1,:,:,:,lb)-ion_face_ym(6,:,:,:,lb)-ion_face_ym(7,:,:,:,lb)-ion_face_ym(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_yp = (ion_face_yp(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_yp(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_yp(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_yp(1,:,:,:,lb)-ion_face_yp(6,:,:,:,lb)-ion_face_yp(7,:,:,:,lb)-ion_face_yp(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_zm = (ion_face_zm(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_zm(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_zm(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_zm(1,:,:,:,lb)-ion_face_zm(6,:,:,:,lb)-ion_face_zm(7,:,:,:,lb)-ion_face_zm(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_zp = (ion_face_zp(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_zp(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_zp(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_zp(1,:,:,:,lb)-ion_face_zp(6,:,:,:,lb)-ion_face_zp(7,:,:,:,lb)-ion_face_zp(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_ei = (ion_face_ei(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_ei(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_ei(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_ei(1,:,:,:,lb)-ion_face_ei(6,:,:,:,lb)-ion_face_ei(7,:,:,:,lb)-ion_face_ei(8,:,:,:,lb))*rho_sc/(M_CO2*amu))
          ion_nd_face_eo = (ion_face_eo(6,:,:,:,lb)*rho_sc/(M_H*amu)) + (ion_face_eo(7,:,:,:,lb)*rho_sc/(M_O*amu)) + (ion_face_eo(8,:,:,:,lb)*rho_sc/(M_O2*amu)) + ((ion_face_eo(1,:,:,:,lb)-ion_face_eo(6,:,:,:,lb)-ion_face_eo(7,:,:,:,lb)-ion_face_eo(8,:,:,:,lb))*rho_sc/(M_CO2*amu))

           ele_face_xm(2,:,:,:,lb) = ion_face_xm(2,:,:,:,lb)-((J_face_xm(1,:,:,:,lb))/(ion_nd_face_xm*q*a_sc))
           ele_face_xp(2,:,:,:,lb) = ion_face_xp(2,:,:,:,lb)-((J_face_xp(1,:,:,:,lb))/(ion_nd_face_xp*q*a_sc))
           ele_face_ym(2,:,:,:,lb) = ion_face_ym(2,:,:,:,lb)-((J_face_ym(1,:,:,:,lb))/(ion_nd_face_ym*q*a_sc))
           ele_face_yp(2,:,:,:,lb) = ion_face_yp(2,:,:,:,lb)-((J_face_yp(1,:,:,:,lb))/(ion_nd_face_yp*q*a_sc))
           ele_face_zm(2,:,:,:,lb) = ion_face_zm(2,:,:,:,lb)-((J_face_zm(1,:,:,:,lb))/(ion_nd_face_zm*q*a_sc))
           ele_face_zp(2,:,:,:,lb) = ion_face_zp(2,:,:,:,lb)-((J_face_zp(1,:,:,:,lb))/(ion_nd_face_zp*q*a_sc))
           ele_face_ei(2,:,:,:,lb) = ion_face_ei(2,:,:,:,lb)-((J_face_ei(1,:,:,:,lb))/(ion_nd_face_ei*q*a_sc))
           ele_face_eo(2,:,:,:,lb) = ion_face_eo(2,:,:,:,lb)-((J_face_eo(1,:,:,:,lb))/(ion_nd_face_eo*q*a_sc))
           ele_face_xm(3,:,:,:,lb) = ion_face_xm(3,:,:,:,lb)-((J_face_xm(2,:,:,:,lb))/(ion_nd_face_xm*q*a_sc))
           ele_face_xp(3,:,:,:,lb) = ion_face_xp(3,:,:,:,lb)-((J_face_xp(2,:,:,:,lb))/(ion_nd_face_xp*q*a_sc))
           ele_face_ym(3,:,:,:,lb) = ion_face_ym(3,:,:,:,lb)-((J_face_ym(2,:,:,:,lb))/(ion_nd_face_ym*q*a_sc))
           ele_face_yp(3,:,:,:,lb) = ion_face_yp(3,:,:,:,lb)-((J_face_yp(2,:,:,:,lb))/(ion_nd_face_yp*q*a_sc))
           ele_face_zm(3,:,:,:,lb) = ion_face_zm(3,:,:,:,lb)-((J_face_zm(2,:,:,:,lb))/(ion_nd_face_zm*q*a_sc))
           ele_face_zp(3,:,:,:,lb) = ion_face_zp(3,:,:,:,lb)-((J_face_zp(2,:,:,:,lb))/(ion_nd_face_zp*q*a_sc))
           ele_face_ei(3,:,:,:,lb) = ion_face_ei(3,:,:,:,lb)-((J_face_ei(2,:,:,:,lb))/(ion_nd_face_ei*q*a_sc))
           ele_face_eo(3,:,:,:,lb) = ion_face_eo(3,:,:,:,lb)-((J_face_eo(2,:,:,:,lb))/(ion_nd_face_eo*q*a_sc))
           ele_face_xm(4,:,:,:,lb) = ion_face_xm(4,:,:,:,lb)-((J_face_xm(3,:,:,:,lb))/(ion_nd_face_xm*q*a_sc))
           ele_face_xp(4,:,:,:,lb) = ion_face_xp(4,:,:,:,lb)-((J_face_xp(3,:,:,:,lb))/(ion_nd_face_xp*q*a_sc))
           ele_face_ym(4,:,:,:,lb) = ion_face_ym(4,:,:,:,lb)-((J_face_ym(3,:,:,:,lb))/(ion_nd_face_ym*q*a_sc))
           ele_face_yp(4,:,:,:,lb) = ion_face_yp(4,:,:,:,lb)-((J_face_yp(3,:,:,:,lb))/(ion_nd_face_yp*q*a_sc))
           ele_face_zm(4,:,:,:,lb) = ion_face_zm(4,:,:,:,lb)-((J_face_zm(3,:,:,:,lb))/(ion_nd_face_zm*q*a_sc))
           ele_face_zp(4,:,:,:,lb) = ion_face_zp(4,:,:,:,lb)-((J_face_zp(3,:,:,:,lb))/(ion_nd_face_zp*q*a_sc))
           ele_face_ei(4,:,:,:,lb) = ion_face_ei(4,:,:,:,lb)-((J_face_ei(3,:,:,:,lb))/(ion_nd_face_ei*q*a_sc))
           ele_face_eo(4,:,:,:,lb) = ion_face_eo(4,:,:,:,lb)-((J_face_eo(3,:,:,:,lb))/(ion_nd_face_eo*q*a_sc))

        else

           ele_face_xm(2,:,:,:,lb) = ion_face_xm(2,:,:,:,lb)
           ele_face_xp(2,:,:,:,lb) = ion_face_xp(2,:,:,:,lb)
           ele_face_ym(2,:,:,:,lb) = ion_face_ym(2,:,:,:,lb)
           ele_face_yp(2,:,:,:,lb) = ion_face_yp(2,:,:,:,lb)
           ele_face_zm(2,:,:,:,lb) = ion_face_zm(2,:,:,:,lb)
           ele_face_zp(2,:,:,:,lb) = ion_face_zp(2,:,:,:,lb)
           ele_face_ei(2,:,:,:,lb) = ion_face_ei(2,:,:,:,lb)
           ele_face_eo(2,:,:,:,lb) = ion_face_eo(2,:,:,:,lb)
           ele_face_xm(3,:,:,:,lb) = ion_face_xm(3,:,:,:,lb)
           ele_face_xp(3,:,:,:,lb) = ion_face_xp(3,:,:,:,lb)
           ele_face_ym(3,:,:,:,lb) = ion_face_ym(3,:,:,:,lb)
           ele_face_yp(3,:,:,:,lb) = ion_face_yp(3,:,:,:,lb)
           ele_face_zm(3,:,:,:,lb) = ion_face_zm(3,:,:,:,lb)
           ele_face_zp(3,:,:,:,lb) = ion_face_zp(3,:,:,:,lb)
           ele_face_ei(3,:,:,:,lb) = ion_face_ei(3,:,:,:,lb)
           ele_face_eo(3,:,:,:,lb) = ion_face_eo(3,:,:,:,lb)
           ele_face_xm(4,:,:,:,lb) = ion_face_xm(4,:,:,:,lb)
           ele_face_xp(4,:,:,:,lb) = ion_face_xp(4,:,:,:,lb)
           ele_face_ym(4,:,:,:,lb) = ion_face_ym(4,:,:,:,lb)
           ele_face_yp(4,:,:,:,lb) = ion_face_yp(4,:,:,:,lb)
           ele_face_zm(4,:,:,:,lb) = ion_face_zm(4,:,:,:,lb)
           ele_face_zp(4,:,:,:,lb) = ion_face_zp(4,:,:,:,lb)
           ele_face_ei(4,:,:,:,lb) = ion_face_ei(4,:,:,:,lb)
           ele_face_eo(4,:,:,:,lb) = ion_face_eo(4,:,:,:,lb)

        end if

! Convert primitive variables to conservative variables -----------------------!                                                           
        call  ele_primitive_to_conserv(lb, ele_face_xm(:,:,:,:,lb), ele_cons_face_xm(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_xp(:,:,:,:,lb), ele_cons_face_xp(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_ym(:,:,:,:,lb), ele_cons_face_ym(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_yp(:,:,:,:,lb), ele_cons_face_yp(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_zm(:,:,:,:,lb), ele_cons_face_zm(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_zp(:,:,:,:,lb), ele_cons_face_zp(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_ei(:,:,:,:,lb), ele_cons_face_ei(:,:,:,:,lb))
        call  ele_primitive_to_conserv(lb, ele_face_eo(:,:,:,:,lb), ele_cons_face_eo(:,:,:,:,lb))

! Compute the average face values ---------------------------------------------!
        ele_face_xm_av(:,:,:,:,lb)=2.
        ele_face_xp_av(:,:,:,:,lb)=2.
        ele_face_ym_av(:,:,:,:,lb)=2.
        ele_face_yp_av(:,:,:,:,lb)=2.
        ele_face_zm_av(:,:,:,:,lb)=2.
        ele_face_zp_av(:,:,:,:,lb)=2.
        ele_face_ex_av(:,:,:,:,lb)=2.

        do i=nguard+1,nxb+nguard
           do j=nguard+1,nyb+nguard
              do k=nguard+1,nzb+nguard
                 ele_face_xm_av(:,i,j,k,lb)=0.5*(ele_face_xm(:,i,j,k,lb)+ele_face_xp(:,i-1,j,k,lb))
                 ele_face_xp_av(:,i,j,k,lb)=0.5*(ele_face_xp(:,i,j,k,lb)+ele_face_xm(:,i+1,j,k,lb))
                 ele_face_ym_av(:,i,j,k,lb)=0.5*(ele_face_ym(:,i,j,k,lb)+ele_face_yp(:,i,j-1,k,lb))
                 ele_face_yp_av(:,i,j,k,lb)=0.5*(ele_face_yp(:,i,j,k,lb)+ele_face_ym(:,i,j+1,k,lb))
                 ele_face_zm_av(:,i,j,k,lb)=0.5*(ele_face_zm(:,i,j,k,lb)+ele_face_zp(:,i,j,k-1,lb))
                 ele_face_zp_av(:,i,j,k,lb)=0.5*(ele_face_zp(:,i,j,k,lb)+ele_face_zm(:,i,j,k+1,lb))
                 ele_face_ex_av(:,i,j,k,lb)=0.5*(ele_face_ei(:,i,j,k,lb)+ele_face_eo(:,i,j,k,lb))

              end do
           end do
        end do

        end subroutine ele_face_surface

!------------------------------------------------------------------------!   
! subroutine ele_residual_normal                                         !
!------------------------------------------------------------------------! 

        subroutine ele_residual_normal(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k, ii, jj, kk
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xm, flx_ym, flx_zm
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xp, flx_yp, flx_zp 
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_tx, flx_ty, flx_tz
        real, dimension(1:ele_nvar,1:nxb,1:nyb,1:nzb) :: flx1, flx2, flx3, flx4, flx5, flx6
        real :: dx, lambda_x, lambda_y, lambda_z, dt_cell, dt_v, ele_p, Vx1, Vx2, Vy3, Vy4, Vz5, Vz6, pv_term, grad_V

! initialize variables --------------------------------------------------!
        dx = (bsize(1,lb)*length_unit/(l_sc*nxb))

! compute face fluxes ---------------------------------------------------!
        call ele_flux_function(ele_face_xm(:,:,:,:,lb), flx_xm, flx_ty, flx_tz)
        call ele_flux_function(ele_face_xp(:,:,:,:,lb), flx_xp, flx_ty, flx_tz)
        call ele_flux_function(ele_face_ym(:,:,:,:,lb), flx_tx, flx_ym, flx_tz)
        call ele_flux_function(ele_face_yp(:,:,:,:,lb), flx_tx, flx_yp, flx_tz)
        call ele_flux_function(ele_face_zm(:,:,:,:,lb), flx_tx, flx_ty, flx_zm)
        call ele_flux_function(ele_face_zp(:,:,:,:,lb), flx_tx, flx_ty, flx_zp)

! Compute the total flux -------------------------------------------------------!
        do k=1,nzb                                                              ! loop over block's cells
           kk = k+nguard
           do j=1,nyb
              jj = j+nguard
              do i=1,nxb
                 ii= i+nguard
                                                                                ! compute average eignvalues
                 lambda_x = max(ele_lambda_xp_av(ii,jj,kk,lb),ele_lambda_xm_av(ii,jj,kk,lb))
                 lambda_y = max(ele_lambda_yp_av(ii,jj,kk,lb),ele_lambda_ym_av(ii,jj,kk,lb))
                 lambda_z = max(ele_lambda_zp_av(ii,jj,kk,lb),ele_lambda_zm_av(ii,jj,kk,lb))

                 flx1(5,i,j,k)=0.5*(flx_xp(5,ii-1,jj,kk)+flx_xm(5,ii,jj,kk))- & ! compute flux on X-
                               ele_art*lambda_x*(ele_cons_face_xm(5,ii,jj,kk,lb)-ele_cons_face_xp(5,ii-1,jj,kk,lb))
                 flx2(5,i,j,k)=0.5*(flx_xp(5,ii,jj,kk)+flx_xm(5,ii+1,jj,kk))- & ! compute flux on X+
                               ele_art*lambda_x*(ele_cons_face_xm(5,ii+1,jj,kk,lb)-ele_cons_face_xp(5,ii,jj,kk,lb))
                 flx3(5,i,j,k)=0.5*(flx_yp(5,ii,jj-1,kk)+flx_ym(5,ii,jj,kk))- & ! compute flux on Y-
                               ele_art*lambda_y*(ele_cons_face_ym(5,ii,jj,kk,lb)-ele_cons_face_yp(5,ii,jj-1,kk,lb))
                 flx4(5,i,j,k)=0.5*(flx_yp(5,ii,jj,kk)+flx_ym(5,ii,jj+1,kk))- & ! compute flux on Y+
                               ele_art*lambda_y*(ele_cons_face_ym(5,ii,jj+1,kk,lb)-ele_cons_face_yp(5,ii,jj,kk,lb))
                 flx5(5,i,j,k)=0.5*(flx_zp(5,ii,jj,kk-1)+flx_zm(5,ii,jj,kk))- & ! compute flux on Z-
                               ele_art*lambda_z*(ele_cons_face_zm(5,ii,jj,kk,lb)-ele_cons_face_zp(5,ii,jj,kk-1,lb))
                 flx6(5,i,j,k)=0.5*(flx_zp(5,ii,jj,kk)+flx_zm(5,ii,jj,kk+1))- & ! compute flux on Z+
                               ele_art*lambda_z*(ele_cons_face_zm(5,ii,jj,kk+1,lb)-ele_cons_face_zp(5,ii,jj,kk,lb))
                                                                           
! Compute the velocity field gradient ------------------------------------------!
                 ele_p  = G0_ele_prim(5,ii,jj,kk,lb)   
                 Vx1 = ele_face_xm_av(2,ii,jj,kk,lb)
                 Vx2 = ele_face_xp_av(2,ii,jj,kk,lb)
                 Vy3 = ele_face_ym_av(3,ii,jj,kk,lb)
                 Vy4 = ele_face_yp_av(3,ii,jj,kk,lb)
                 Vz5 = ele_face_zm_av(4,ii,jj,kk,lb)
                 Vz6 = ele_face_zp_av(4,ii,jj,kk,lb)
                 grad_V = (Vx2-Vx1+Vy4-Vy3+Vz6-Vz5)/dx
                 pv_term = (gamma-1)*ele_p*grad_V

! Compute electron residual ----------------------------------------------------!
                 ele_R(5,i,j,k,lb) = ele_R(5,i,j,k,lb)+(flx6(5,i,j,k)-flx5(5,i,j,k)+flx4(5,i,j,k)-flx3(5,i,j,k) &
                                                       +flx2(5,i,j,k)-flx1(5,i,j,k))/dx-ele_src(5,i,j,k)+pv_term

                 if (init_sstep) then

                    dt_v = CFL/((2.*(lambda_x+lambda_y+lambda_z)/dx) +(1./cell_time_stif(2,i,j,k,lb))+1.E-50)
                    if (ele_R(5,i,j,k,lb)>0.) dt_v = min(dt_v,fac*(G0_ele_cons(5,ii,jj,kk,lb))/(ele_R(5,i,j,k,lb)+1.E-50)) 
                    
                    dt_cell = min(dt_v,maxtimestep)
                    
                    cell_time(2,i,j,k,lb) = dt_cell                             ! store cell's time step
                    local_time(2,lb) = min(dt_cell, local_time(2,lb))           ! compute block's electron time step 
                    local_time(4,lb) = min(dt_cell, local_time(4,lb))           ! compute block's absolute time step 
 
                 end if
 
              end do
           end do
        end do

        return

        end subroutine ele_residual_normal

!-------------------------------------------------------------------------------!   
! subroutine ele_adv_normal                                                     !
!-------------------------------------------------------------------------------! 

        subroutine ele_adv_normal(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k, ii, jj, kk
        real :: dt
        real, dimension(1:nxb,1:nyb,1:nzb) :: res

! update the solution ----------------------------------------------------------!
        do k = nguard+1,nzb+nguard
           kk = k-nguard
           do j = nguard+1,nyb+nguard
              jj = j-nguard
              do i = nguard+1,nxb+nguard 
                 ii = i-nguard 
                 
                 if ((mixed_time_stepping).and.(global_time_step.lt.istep)) then
                    dt = dtsubstep                                              ! use global time stepping
                 else
                    dt = cell_subtime(ii,jj,kk,lb)                              ! use local time stepping 
                 end if

                 G1_ele_cons(5,i,j,k,lb)=G0_ele_cons(5,i,j,k,lb)-(dt*ele_R(5,ii,jj,kk,lb))

                 res(ii,jj,kk)=abs(ele_R(5,ii,jj,kk,lb))  

              end do
           end do
        end do         

! compute the local residual ---------------------------------------------------!
        local_residu(2,lb) = maxval(res)

        return

        end subroutine ele_adv_normal

!-------------------------------------------------------------------------------!   
! subroutine ele_residual_surface                                               !
!-------------------------------------------------------------------------------! 

        subroutine ele_residual_surface(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep

        integer :: i, j, k, ii, jj, kk, face
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xm, flx_ym, flx_zm
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xp, flx_yp, flx_zp 
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_tx, flx_ty, flx_tz
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_eo_xm, flx_eo_ym, flx_eo_zm
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_ei_xp, flx_ei_yp, flx_ei_zp
        real, dimension(1:ele_nvar,1:nxb,1:nyb,1:nzb) :: flx1, flx2, flx3, flx4, flx5, flx6, flx7, flxa, flxb, flxc, flx, flx_tot
        real :: charac, lambda_x_e, lambda_y_e, lambda_z_e, lambda_x, lambda_y, lambda_z
        real :: ele_p, surface, nx, ny, nz, volume, Vx1, Vx2, Vy3, Vy4, Vz5, Vz6, pv_term, grad_V
        real :: dt_cell, dt_v, dx_scale, dx


 
! Compute cell size -------- ---------------------------------------------------!
        dx = (bsize(1,lb)*length_unit/(l_sc*nxb))
        dx_scale = length_unit/l_sc

! Compute face fluxes ----------------------------------------------------------!
        call ele_flux_function(ele_face_xm(:,:,:,:,lb), flx_xm, flx_ty, flx_tz)
        call ele_flux_function(ele_face_xp(:,:,:,:,lb), flx_xp, flx_ty, flx_tz)
        call ele_flux_function(ele_face_ym(:,:,:,:,lb), flx_tx, flx_ym, flx_tz)
        call ele_flux_function(ele_face_yp(:,:,:,:,lb), flx_tx, flx_yp, flx_tz)
        call ele_flux_function(ele_face_zm(:,:,:,:,lb), flx_tx, flx_ty, flx_zm)
        call ele_flux_function(ele_face_zp(:,:,:,:,lb), flx_tx, flx_ty, flx_zp)
        call ele_flux_function(ele_face_eo(:,:,:,:,lb), flx_eo_xm, flx_eo_ym, flx_eo_zm)
        call ele_flux_function(ele_face_ei(:,:,:,:,lb), flx_ei_xp, flx_ei_yp, flx_ei_zp)

! Compute the flux functions ---------------------------------------------------!
        do k=1,nzb                                                              ! loop over block's cells
           kk = k+nguard
           do j=1,nyb
              jj = j+nguard
              do i=1,nxb
                 ii= i+nguard

! Case of a non surface cell ---------------------------------------------------!
                 if (loc_cell_flags(ii,jj,kk).eq.0) then
                                                                                ! compute average eignvalues
                    lambda_x = max(ele_lambda_xp_av(ii,jj,kk,lb),ele_lambda_xm_av(ii,jj,kk,lb))
                    lambda_y = max(ele_lambda_yp_av(ii,jj,kk,lb),ele_lambda_ym_av(ii,jj,kk,lb))
                    lambda_z = max(ele_lambda_zp_av(ii,jj,kk,lb),ele_lambda_zm_av(ii,jj,kk,lb))
                    
                    flx1(5,i,j,k)=0.5*(flx_xp(5,ii-1,jj,kk)+flx_xm(5,ii,jj,kk))-& ! compute flux on X-
                         ele_art*lambda_x*(ele_cons_face_xm(5,ii,jj,kk,lb)-ele_cons_face_xp(5,ii-1,jj,kk,lb))
                    flx2(5,i,j,k)=0.5*(flx_xp(5,ii,jj,kk)+flx_xm(5,ii+1,jj,kk))-& ! compute flux on X+
                         ele_art*lambda_x*(ele_cons_face_xm(5,ii+1,jj,kk,lb)-ele_cons_face_xp(5,ii,jj,kk,lb))
                    flx3(5,i,j,k)=0.5*(flx_yp(5,ii,jj-1,kk)+flx_ym(5,ii,jj,kk))-& ! compute flux on Y-
                         ele_art*lambda_y*(ele_cons_face_ym(5,ii,jj,kk,lb)-ele_cons_face_yp(5,ii,jj-1,kk,lb))
                    flx4(5,i,j,k)=0.5*(flx_yp(5,ii,jj,kk)+flx_ym(5,ii,jj+1,kk))-& ! compute flux on Y+
                         ele_art*lambda_y*(ele_cons_face_ym(5,ii,jj+1,kk,lb)-ele_cons_face_yp(5,ii,jj,kk,lb))
                    flx5(5,i,j,k)=0.5*(flx_zp(5,ii,jj,kk-1)+flx_zm(5,ii,jj,kk))-& ! compute flux on Z-
                         ele_art*lambda_z*(ele_cons_face_zm(5,ii,jj,kk,lb)-ele_cons_face_zp(5,ii,jj,kk-1,lb))
                    flx6(5,i,j,k)=0.5*(flx_zp(5,ii,jj,kk)+flx_zm(5,ii,jj,kk+1))-& ! compute flux on Z+
                         ele_art*lambda_z*(ele_cons_face_zm(5,ii,jj,kk+1,lb)-ele_cons_face_zp(5,ii,jj,kk,lb)) 

! Compute the velocity field gradient ------------------------------------------!
                    ele_p  = G0_ele_prim(5,ii,jj,kk,lb)   
                    Vx1 = ele_face_xm_av(2,ii,jj,kk,lb)
                    Vx2 = ele_face_xp_av(2,ii,jj,kk,lb)
                    Vy3 = ele_face_ym_av(3,ii,jj,kk,lb)
                    Vy4 = ele_face_yp_av(3,ii,jj,kk,lb)
                    Vz5 = ele_face_zm_av(4,ii,jj,kk,lb)
                    Vz6 = ele_face_zp_av(4,ii,jj,kk,lb)
                    grad_V = (Vx2-Vx1+Vy4-Vy3+Vz6-Vz5)/dx
                    pv_term = (gamma-1)*ele_p*grad_V
                    
! Compute electron residual ----------------------------------------------------!
                    ele_R(5,i,j,k,lb) = ele_R(5,i,j,k,lb)+(flx6(5,i,j,k)-flx5(5,i,j,k)+flx4(5,i,j,k)-flx3(5,i,j,k) &
                                                         +flx2(5,i,j,k)-flx1(5,i,j,k))/dx-ele_src(5,i,j,k)+pv_term

                    if (init_sstep) then

                       dt_v = CFL/((2.*(lambda_x+lambda_y+lambda_z)/dx)+(1./cell_time_stif(2,i,j,k,lb))+1.E-50)
                       if (ele_R(5,i,j,k,lb)>0.) dt_v = min(dt_v,fac*(G0_ele_cons(5,ii,jj,kk,lb))/(ele_R(5,i,j,k,lb)+1.E-50)) 
                       
                       dt_cell = min(dt_v,maxtimestep)
                    
                       cell_time(2,i,j,k,lb) = dt_cell                                ! store cell's time step
                       local_time(2,lb) = min(dt_cell, local_time(2,lb))              ! compute block's electron time step 
                       local_time(4,lb) = min(dt_cell, local_time(4,lb))              ! compute block's absolute time step  

                    end if
                 end if

! Case of a surface cell -------------------------------------------------------!
                 if (loc_cell_flags(ii,jj,kk).gt.0) then
                   
                                                                                ! compute average eignvalues 
                    lambda_x = 0.
                    lambda_y = 0.
                    lambda_z = 0.
                    if (surf_face_data(3,1,ii,jj,kk,bflags(1,lb)).ne.0) lambda_x = max(lambda_x,ele_lambda_xm_av(ii,jj,kk,lb))
                    if (surf_face_data(3,2,ii,jj,kk,bflags(1,lb)).ne.0) lambda_x = max(lambda_x,ele_lambda_xp_av(ii,jj,kk,lb))
                    if (surf_face_data(3,3,ii,jj,kk,bflags(1,lb)).ne.0) lambda_y = max(lambda_y,ele_lambda_ym_av(ii,jj,kk,lb))
                    if (surf_face_data(3,4,ii,jj,kk,bflags(1,lb)).ne.0) lambda_y = max(lambda_y,ele_lambda_yp_av(ii,jj,kk,lb))
                    if (surf_face_data(3,5,ii,jj,kk,bflags(1,lb)).ne.0) lambda_z = max(lambda_z,ele_lambda_zm_av(ii,jj,kk,lb))
                    if (surf_face_data(3,6,ii,jj,kk,bflags(1,lb)).ne.0) lambda_z = max(lambda_z,ele_lambda_zp_av(ii,jj,kk,lb))
                    if (surf_face_data(3,7,ii,jj,kk,bflags(1,lb)).ne.0) then
                       lambda_x = max(lambda_x,ele_lambda_ex_av(ii,jj,kk,lb))
                       lambda_y = max(lambda_y,ele_lambda_ey_av(ii,jj,kk,lb))
                       lambda_z = max(lambda_z,ele_lambda_ez_av(ii,jj,kk,lb))
                    end if

                    charac = 0.
                    flx_tot(:,:,:,:) = 0.
                    grad_V = 0.

                    do face=1,7                                                 ! loop over cell's faces
                       if (surf_face_data(3,face,ii,jj,kk,bflags(1,lb)).ne.0) then
                          select case(face)
                             case(1)                                            ! compute flux on X-
                                flx1(5,i,j,k)=0.5*(flx_xp(5,ii-1,jj,kk)+flx_xm(5,ii,jj,kk))- &
                                     ele_art*lambda_x*(ele_cons_face_xm(5,ii,jj,kk,lb)-ele_cons_face_xp(5,ii-1,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)-(flx1(5,i,j,k)*surface)
                                charac = charac+lambda_x*surface
                                grad_V =grad_V-(ele_face_xm_av(2,ii,jj,kk,lb)*surface)
                             case(2)                                            ! compute flux on X+
                                flx2(5,i,j,k)=0.5*(flx_xp(5,ii,jj,kk)+flx_xm(5,ii+1,jj,kk))- &
                                     ele_art*lambda_x*(ele_cons_face_xm(5,ii+1,jj,kk,lb)-ele_cons_face_xp(5,ii,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)+(flx2(5,i,j,k)*surface)
                                charac = charac+lambda_x*surface
                                grad_V =grad_V+(ele_face_xp_av(2,ii,jj,kk,lb)*surface)
                             case(3)                                            ! compute flux on Y-
                                flx3(5,i,j,k)=0.5*(flx_yp(5,ii,jj-1,kk)+flx_ym(5,ii,jj,kk))- &
                                     ele_art*lambda_y*(ele_cons_face_ym(5,ii,jj,kk,lb)-ele_cons_face_yp(5,ii,jj-1,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)-(flx3(5,i,j,k)*surface)
                                charac = charac+lambda_y*surface
                                grad_V =grad_V-(ele_face_ym_av(3,ii,jj,kk,lb)*surface)
                                      case(4)                                            ! compute flux on Y+
                                flx4(5,i,j,k)=0.5*(flx_yp(5,ii,jj,kk)+flx_ym(5,ii,jj+1,kk))- &
                                     ele_art*lambda_y*(ele_cons_face_ym(5,ii,jj+1,kk,lb)-ele_cons_face_yp(5,ii,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)+(flx4(5,i,j,k)*surface)
                                charac = charac+lambda_y*surface
                                grad_V =grad_V+(ele_face_yp_av(3,ii,jj,kk,lb)*surface)
                             case(5)                                            ! compute flux on Z-
                                flx5(5,i,j,k)=0.5*(flx_zp(5,ii,jj,kk-1)+flx_zm(5,ii,jj,kk))- &
                                     ele_art*lambda_z*(ele_cons_face_zm(5,ii,jj,kk,lb)-ele_cons_face_zp(5,ii,jj,kk-1,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)-(flx5(5,i,j,k)*surface)
                                charac = charac+lambda_z*surface
                                grad_V =grad_V-(ele_face_zm_av(4,ii,jj,kk,lb)*surface)
                             case(6)                                            ! compute flux on Z+ 
                                flx6(5,i,j,k)=0.5*(flx_zp(5,ii,jj,kk)+flx_zm(5,ii,jj,kk+1))- &
                                     ele_art*lambda_z*(ele_cons_face_zm(5,ii,jj,kk+1,lb)-ele_cons_face_zp(5,ii,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)+(flx6(5,i,j,k)*surface)
                                charac = charac+lambda_z*surface
                                grad_V =grad_V+(ele_face_zp_av(4,ii,jj,kk,lb)*surface)
                             case(7)                                            ! compute flux on Cut face
                                nx = surf_face_data(4,face,ii,jj,kk,bflags(1,lb))
                                ny = surf_face_data(5,face,ii,jj,kk,bflags(1,lb))
                                nz = surf_face_data(6,face,ii,jj,kk,bflags(1,lb))
                                flxa(5,i,j,k)=0.5*(flx_ei_xp(5,ii,jj,kk)+flx_eo_xm(5,ii,jj,kk))- &
                                     sign(1.,nx)*ele_art*lambda_x*(ele_cons_face_eo(5,ii,jj,kk,lb)-ele_cons_face_ei(5,ii,jj,kk,lb))
                                flxb(5,i,j,k)=0.5*(flx_ei_yp(5,ii,jj,kk)+flx_eo_ym(5,ii,jj,kk))- &
                                     sign(1.,ny)*ele_art*lambda_y*(ele_cons_face_eo(5,ii,jj,kk,lb)-ele_cons_face_ei(5,ii,jj,kk,lb))
                                flxc(5,i,j,k)=0.5*(flx_ei_zp(5,ii,jj,kk)+flx_eo_zm(5,ii,jj,kk))- &
                                     sign(1.,nz)*ele_art*lambda_z*(ele_cons_face_eo(5,ii,jj,kk,lb)-ele_cons_face_ei(5,ii,jj,kk,lb))
                                flx7(5,i,j,k) = (flxa(5,i,j,k)*nx+flxb(5,i,j,k)*ny+flxc(5,i,j,k)*nz)
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(5,i,j,k)=flx_tot(5,i,j,k)+(flx7(5,i,j,k)*surface)
                                charac = charac+(lambda_x+lambda_y+lambda_z)*surface
                                grad_V =grad_V+(ele_face_ex_av(2,ii,jj,kk,lb)*nx+ele_face_ex_av(3,ii,jj,kk,lb)*ny+ &
                                                ele_face_ex_av(4,ii,jj,kk,lb)*nz)*surface
                          end select
                       end if
                    end do
                                                                                ! compute cell's volume
                    volume = surf_cell_data(1,ii,jj,kk,bflags(1,lb))*(dx_scale**3.)

! Compute the velocity field gradient ------------------------------------------!
                    ele_p  = G0_ele_prim(5,ii,jj,kk,lb)   
                    grad_V=grad_V/volume
                    pv_term = (gamma-1)*ele_p*grad_V

! Compute ele residual ---------------------------------------------------------!
                    ele_R(5,i,j,k,lb)=(flx_tot(5,i,j,k)/volume)-ele_src(5,i,j,k)+pv_term

                    if (init_sstep) then

                       dt_v = (CFL*volume)/(charac+(volume/cell_time_stif(2,i,j,k,lb))+1.e-50)
                       if (ele_R(5,i,j,k,lb)>0.) dt_v = min(dt_v,fac*(G0_ele_cons(5,ii,jj,kk,lb))/(ele_R(5,i,j,k,lb)+1.E-50))

                       dt_cell = min(dt_v,maxtimestep)

                       cell_time(2,i,j,k,lb) = dt_cell                                ! store cell's time step
                       local_time(2,lb) = min(dt_cell, local_time(2,lb))              ! compute block's electron time step 

                    end if
                 end if
          
              end do
           end do
        end do

        return
        
        end subroutine ele_residual_surface

!------------------------------------------------------------------------!
! subroutine ele_adv_surface                                             !
!------------------------------------------------------------------------! 

        subroutine ele_adv_surface(mype, lb, istep)

        implicit none

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k, ii, jj, kk
        real :: dt
        real, dimension(1:nxb,1:nyb,1:nzb) :: res

        res(:,:,:) = 0.
        G1_ele_cons(5,:,:,:,lb) = 2. 

! update the solution ---------------------------------------------------!
        do k = nguard+1,nzb+nguard
           kk = k-nguard
           do j = nguard+1,nyb+nguard
              jj = j-nguard
              do i = nguard+1,nxb+nguard 
                 ii = i-nguard 

                 if (loc_cell_flags(i,j,k).ge.0) then                           ! check if it is not a core cell

                    if (loc_cell_flags(i,j,k).gt.0) then
                       dt = cell_subtime(ii,jj,kk,lb) 
                    else
                       if ((mixed_time_stepping).and.(global_time_step.lt.istep)) then
                          dt = dtsubstep                                        ! use global time stepping
                       else
                          dt = cell_subtime(ii,jj,kk,lb)                        ! use local time stepping
                       end if
                    end if
                                                                                ! compute next solution
                    G1_ele_cons(5,i,j,k,lb)=G0_ele_cons(5,i,j,k,lb)-(dt*ele_R(5,ii,jj,kk,lb))

                    res(ii,jj,kk)=abs(ele_R(5,ii,jj,kk,lb))

                 end if

              end do
           end do
        end do          
                       
! compute the local residual --------------------------------------------!
        local_residu(2,lb) = maxval(res)
        
        return

        end subroutine ele_adv_surface
!------------------------------------------------------------------------!   
! subroutine ele_interpol_to_faces_normal                                !
!------------------------------------------------------------------------! 

        subroutine ele_interpol_to_faces_normal(lb, istep, vec, fxm, fxp, fym, fyp, fzm, fzp, ord)

        implicit none 
        
        integer, intent(in) :: lb, istep, ord
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: vec
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: fxm, fxp, fym, fyp, fzm, fzp
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: fei
        integer :: i, j, k

        fxm(1,:,:,:) = 1.
        fxp(1,:,:,:) = 1.
        fym(1,:,:,:) = 1.
        fyp(1,:,:,:) = 1.
        fzm(1,:,:,:) = 1.
        fzp(1,:,:,:) = 1.

! reconstruct cell edges ------------------------------------------------!
        do k=2,nzb+2*nguard-1
           do j=2,nyb+2*nguard-1
              do i=2,nxb+2*nguard-1
                 call ele_reconstruct_edge(lb, vec(:,i,j,k), vec(:,i-1,j,k), vec(:,i+1,j,k), vec(:,i,j-1,k), &
                   vec(:,i,j+1,k), vec(:,i,j,k-1), vec(:,i,j,k+1), &
                   vec(:,i-1,j-1,k), vec(:,i+1,j-1,k), vec(:,i+1,j+1,k), vec(:,i-1,j+1,k), &
                   vec(:,i-1,j-1,k-1), vec(:,i+1,j-1,k-1), vec(:,i+1,j+1,k-1), vec(:,i-1,j+1,k-1), &
                   vec(:,i-1,j-1,k+1), vec(:,i+1,j-1,k+1), vec(:,i+1,j+1,k+1), vec(:,i-1,j+1,k+1), &
                   vec(:,i-1,j,k-1), vec(:,i+1,j,k-1), vec(:,i,j-1,k-1), vec(:,i,j+1,k-1), &
                   vec(:,i-1,j,k+1), vec(:,i+1,j,k+1), vec(:,i,j-1,k+1), vec(:,i,j+1,k+1), &
                   fxm(:,i,j,k), fxp(:,i,j,k), fym(:,i,j,k), fyp(:,i,j,k), fzm(:,i,j,k), fzp(:,i,j,k), fei(:,i,j,k), ord, dxcell(i,j,k,lb))
              end do
           end do
        end do

        return
        end subroutine ele_interpol_to_faces_normal

!------------------------------------------------------------------------!   
! subroutine ele_interpol_to_faces_surface                               !
!------------------------------------------------------------------------! 

        subroutine ele_interpol_to_faces_surface(lb, istep, vec, fxm, fxp, fym, fyp, fzm, fzp, fei, feo, ord)

        implicit none 
        
        integer, intent(in) :: lb, istep, ord
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: vec
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: fxm, fxp, fym, fyp, fzm, fzp, feo, fei
        integer :: i, j, k, ii, jj, kk, face, nuorder
        real :: xc, yc, zc

        fxm(1,:,:,:) = 1.
        fxp(1,:,:,:) = 1.
        fym(1,:,:,:) = 1.
        fyp(1,:,:,:) = 1.
        fzm(1,:,:,:) = 1.
        fzp(1,:,:,:) = 1.
        fei(1,:,:,:) = 1.
        feo(1,:,:,:) = 1.
         
! reconstruct cell edges ------------------------------------------------!
        do k=2,nzb+2*nguard-1
           do j=2,nyb+2*nguard-1
              do i=2,nxb+2*nguard-1


                 if (loc_cell_flags(i,j,k).ne.-1) then                          ! check that it is not a nucleus cell

                    if(loc_cell_flags(i,j,k).eq.0) then                         ! case of non surface cell           

                       nuorder = ord     
                       call ele_reconstruct_edge(lb, vec(:,i,j,k), vec(:,i-1,j,k), vec(:,i+1,j,k), vec(:,i,j-1,k), &
                       vec(:,i,j+1,k), vec(:,i,j,k-1), vec(:,i,j,k+1), &
                       vec(:,i-1,j-1,k), vec(:,i+1,j-1,k), vec(:,i+1,j+1,k), vec(:,i-1,j+1,k), &
                       vec(:,i-1,j-1,k-1), vec(:,i+1,j-1,k-1), vec(:,i+1,j+1,k-1), vec(:,i-1,j+1,k-1), &
                       vec(:,i-1,j-1,k+1), vec(:,i+1,j-1,k+1), vec(:,i+1,j+1,k+1), vec(:,i-1,j+1,k+1), &
                       vec(:,i-1,j,k-1), vec(:,i+1,j,k-1), vec(:,i,j-1,k-1), vec(:,i,j+1,k-1), &
                       vec(:,i-1,j,k+1), vec(:,i+1,j,k+1), vec(:,i,j-1,k+1), vec(:,i,j+1,k+1), &
                       fxm(:,i,j,k), fxp(:,i,j,k), fym(:,i,j,k), fyp(:,i,j,k), fzm(:,i,j,k), fzp(:,i,j,k), fei(:,i,j,k), nuorder, dxcell(i,j,k,lb))

                    end if

                    if (loc_cell_flags(i,j,k).gt.0) then                        ! case of a surface cell

                       nuorder = 2   
                       call ele_reconstruct_edge2(lb, i, j, k, vec(:,i,j,k), vec(:,i-1,j,k), vec(:,i+1,j,k), vec(:,i,j-1,k), &
                       vec(:,i,j+1,k), vec(:,i,j,k-1), vec(:,i,j,k+1), &
                       vec(:,i-1,j-1,k), vec(:,i+1,j-1,k), vec(:,i+1,j+1,k), vec(:,i-1,j+1,k), &
                       vec(:,i-1,j-1,k-1), vec(:,i+1,j-1,k-1), vec(:,i+1,j+1,k-1), vec(:,i-1,j+1,k-1), &
                       vec(:,i-1,j-1,k+1), vec(:,i+1,j-1,k+1), vec(:,i+1,j+1,k+1), vec(:,i-1,j+1,k+1), &
                       vec(:,i-1,j,k-1), vec(:,i+1,j,k-1), vec(:,i,j-1,k-1), vec(:,i,j+1,k-1), &
                       vec(:,i-1,j,k+1), vec(:,i+1,j,k+1), vec(:,i,j-1,k+1), vec(:,i,j+1,k+1), &
                       fxm(:,i,j,k), fxp(:,i,j,k), fym(:,i,j,k), fyp(:,i,j,k), fzm(:,i,j,k), fzp(:,i,j,k), fei(:,i,j,k),  nuorder)

                       if (surf_face_data(3,7,i,j,k,bflags(1,lb)).ne.0) then    ! apply boundary condition
                          xc = surf_face_data(8,7,i,j,k,bflags(1,lb))*length_unit
                          yc = surf_face_data(9,7,i,j,k,bflags(1,lb))*length_unit
                          zc = surf_face_data(10,7,i,j,k,bflags(1,lb))*length_unit
                          call ele_set_surf_bc(lb, xc, yc, zc, fei(:,i,j,k), feo(:,i,j,k))
                       end if

                       if (surf_face_data(3,1,i,j,k,bflags(1,lb)).eq.0) then    ! make sure that all arrays have valid values
                          fxp(:,i-1,j,k) = vec(:,i,j,k)
                          fxm(:,i,j,k) = vec(:,i,j,k)
                       end if
                       if (surf_face_data(3,2,i,j,k,bflags(1,lb)).eq.0) then
                          fxm(:,i+1,j,k) = vec(:,i,j,k)
                          fxp(:,i,j,k) = vec(:,i,j,k)
                       end if
                       if (surf_face_data(3,3,i,j,k,bflags(1,lb)).eq.0) then
                          fyp(:,i,j-1,k) = vec(:,i,j,k)
                          fym(:,i,j,k) = vec(:,i,j,k)
                       end if
                       if (surf_face_data(3,4,i,j,k,bflags(1,lb)).eq.0) then
                          fym(:,i,j+1,k) = vec(:,i,j,k)
                          fyp(:,i,j,k) = vec(:,i,j,k)
                       end if
                       if (surf_face_data(3,5,i,j,k,bflags(1,lb)).eq.0) then
                          fzp(:,i,j,k-1) = vec(:,i,j,k)
                          fzm(:,i,j,k) = vec(:,i,j,k)
                       end if
                       if (surf_face_data(3,6,i,j,k,bflags(1,lb)).eq.0) then
                          fzm(:,i,j,k+1) = vec(:,i,j,k)
                          fzp(:,i,j,k) = vec(:,i,j,k)
                       end if
                       if (surf_face_data(3,7,i,j,k,bflags(1,lb)).eq.0) then
                          feo(:,i,j,k) = vec(:,i,j,k)
                          fei(:,i,j,k) = vec(:,i,j,k)
                       end if
                    end if
                       
                 end if
              end do
           end do
        end do 
        
        return
        end subroutine  ele_interpol_to_faces_surface

!------------------------------------------------------------------------!   
! subroutine ele_set_surf_bc                                             !
!------------------------------------------------------------------------!  
        subroutine ele_set_surf_bc(lb, xc, yc, zc, fi, fo)

        implicit none 
        integer, intent(in) :: lb
        real, intent(in) :: xc, yc, zc
        real, dimension(1:ele_nvar), intent(inout) :: fi, fo
        real ::  nr, nt, np, r, theta, phi, vx, vy, vz

        r = sqrt(xc**2.+yc**2.+zc**2.)
        theta = acos(zc/r)
        phi = atan2(yc,xc)

        fo(:) = fi(:)

        vx = fi(2)
        vy = fi(3)
        vz = fi(4)

        nr = (sin(theta)*cos(phi)*vx)+(sin(theta)*sin(phi)*vy)+(cos(theta)*vz)
        nt = (cos(theta)*cos(phi)*vx)+(cos(theta)*sin(phi)*vy)-(sin(theta)*vz)
        np = -(sin(phi)*vx)+(cos(phi)*vy)
        if (nr>0.) then
           nr = -nr
        end if
        vx = (sin(theta)*cos(phi)*nr)+(cos(theta)*cos(phi)*nt)-(sin(phi)*np)
        vy = (sin(theta)*sin(phi)*nr)+(cos(theta)*sin(phi)*nt)+(cos(phi)*np)
        vz = (cos(theta)*nr)-(sin(theta)*nt)
 
        fo(2) = vx
        fo(3) = vy
        fo(4) = vz

        return
        end subroutine ele_set_surf_bc

!-------------------------------------------------------------------------------!   
! Subroutine ele_primitive_to_conserv                                           !
!-------------------------------------------------------------------------------!

        subroutine ele_primitive_to_conserv(lb, prim, conserv)

        implicit none
      
        integer, intent(in) :: lb
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in)  :: prim
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out) :: conserv

        conserv(1,:,:,:) = prim(1,:,:,:)
        conserv(2,:,:,:) = prim(1,:,:,:)*prim(2,:,:,:)
        conserv(3,:,:,:) = prim(1,:,:,:)*prim(3,:,:,:)
        conserv(4,:,:,:) = prim(1,:,:,:)*prim(4,:,:,:)
        conserv(5,:,:,:) = prim(5,:,:,:)

        return
        end subroutine ele_primitive_to_conserv

!-------------------------------------------------------------------------------!   
! subroutine ele_conserv_to_primitive                                           !
!-------------------------------------------------------------------------------!

        subroutine ele_conserv_to_primitive(lb, conserv, prim)

        implicit none

        integer, intent(in) :: lb
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: prim
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: conserv

        prim(1,:,:,:) = conserv(1,:,:,:)
        prim(2,:,:,:) = conserv(2,:,:,:)/prim(1,:,:,:)
        prim(3,:,:,:) = conserv(3,:,:,:)/prim(1,:,:,:)
        prim(4,:,:,:) = conserv(4,:,:,:)/prim(1,:,:,:)
        prim(5,:,:,:) = conserv(5,:,:,:)

        return
        end subroutine ele_conserv_to_primitive

!------------------------------------------------------------------------!   
! subroutine ele_adjust_storage      !forcing neutrality                 !
!------------------------------------------------------------------------!

        subroutine ele_adjust_storage(lb, conserv, ion_cons)

        implicit none

        integer, intent(in) :: lb
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(inout) :: conserv
        real, dimension(1:ion_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: ion_cons
        integer :: i ,j, k
	    real :: ele_T

        conserv(1,:,:,:) = (ion_cons(6,:,:,:)*rho_sc)/(M_H*amu*nd_sc) + (ion_cons(7,:,:,:)*rho_sc)/(M_O*amu*nd_sc) + (ion_cons(8,:,:,:)*rho_sc)/(M_O2*amu*nd_sc) + ((ion_cons(1,:,:,:)-ion_cons(6,:,:,:)-ion_cons(7,:,:,:)-ion_cons(8,:,:,:))*rho_sc)/(M_CO2*amu*nd_sc)
        
        do k=nguard+1,nzb+nguard
           do j=nguard+1,nyb+nguard
              do i=nguard+1,nxb+nguard 

                 if (ieee_is_nan(conserv(1,i,j,k))) then
                    exception_flag = 1
                    print*,'Exception detected Ele_Module :', trans_mype, lb, i,j,k, conserv(1,i,j,k)
                 end if 

                 ele_T = (conserv(5,i,j,k)*p_sc/(Boltz*conserv(1,i,j,k)*nd_sc))
                 if (ele_T<ele_T_min) conserv(5,i,j,k)=conserv(1,i,j,k)*nd_sc*boltz*ele_T_min/p_sc
                 if (ele_T>ele_T_max) conserv(5,i,j,k)=conserv(1,i,j,k)*nd_sc*boltz*ele_T_max/p_sc
              end do
           end do
        end do


        return
        end subroutine ele_adjust_storage


!------------------------------------------------------------------------!   
! subroutine ele_flux_function                                           !
!------------------------------------------------------------------------!

        subroutine ele_flux_function(prim, flx_x, flx_y, flx_z)

        implicit none
      
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: prim
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: flx_x, flx_y, flx_z
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: vx, vy, vz, p

        vx(:,:,:)   = prim(2,:,:,:)
        vy(:,:,:)   = prim(3,:,:,:)
        vz(:,:,:)   = prim(4,:,:,:)
        p(:,:,:)    = prim(5,:,:,:)

        flx_x(1:4,:,:,:) = 0.
        flx_y(1:4,:,:,:) = 0.
        flx_z(1:4,:,:,:) = 0.

        flx_x(5,:,:,:) = p(:,:,:)*vx(:,:,:)
        flx_y(5,:,:,:) = p(:,:,:)*vy(:,:,:)
        flx_z(5,:,:,:) = p(:,:,:)*vz(:,:,:)

        return
        end subroutine ele_flux_function

!-------------------------------------------------------------------------------!   
! Subroutine ele_reconstruct_edge                                               !
!-------------------------------------------------------------------------------!

        subroutine ele_reconstruct_edge(lb, W_0, W_1, W_2, W_3, W_4, W_5, W_6, &
        C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12, E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, &
        W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7, ord, dx)

        implicit none

        integer, intent(in) :: ord, lb
        real, intent(in) :: dx
        real, dimension(ele_nvar), intent(in) :: W_0, W_1, W_2, W_3, W_4, W_5, W_6
        real, dimension(ele_nvar), intent(in) :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12
        real, dimension(ele_nvar), intent(in) :: E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8
        real, dimension(ele_nvar), intent(out):: W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7
        real, dimension(6,ele_nvar) :: tmp_pos
        real, dimension(ele_nvar) :: xdUM, xdUP, ydUM, ydUP, zdUM, zdUP  
        real, dimension(ele_nvar) ::  limiter1, limiter, poly, y, sig, dU2, minpath, maxpath
        integer :: i, face
        real :: trigger , phi, k

        W_0_edge7(:) = W_0(:)
        if (ord .eq. 1) then
           W_0_edge1(:) = W_0(:)
           W_0_edge2(:) = W_0(:)
           W_0_edge3(:) = W_0(:)
           W_0_edge4(:) = W_0(:)
           W_0_edge5(:) = W_0(:)
           W_0_edge6(:) = W_0(:)
           return
        end if

        phi = 1. 
        k = 1./3.
        trigger = (0.25*dx)**3.

        xdUP(:) = 0.125*(0.5*(C_3(:)+C_2(:)-W_4(:)-W_3(:)+E_2(:)+E_6(:)-W_6(:)-W_5(:))+ &
                    0.25*(C_6(:)+C_7(:)+C_10(:)+C_11(:)-E_3(:)-E_4(:)-E_7(:)-E_8(:))+W_2(:)-W_0(:))

        ydUP(:) = 0.125*(0.5*(C_3(:)+C_4(:)-W_2(:)-W_1(:)+E_8(:)+E_4(:)-W_5(:)-W_6(:))+ & 
                    0.25*(C_11(:)+C_12(:)+C_7(:)+C_8(:)-E_6(:)-E_5(:)-E_2(:)-E_1(:))+W_4(:)-W_0(:))

        zdUP(:) = 0.125*(0.5*(E_6(:)+E_5(:)-W_2(:)-W_1(:)+E_8(:)+E_7(:)-W_3(:)-W_4(:))+ & 
                    0.25*(C_11(:)+C_12(:)+C_10(:)+C_9(:)-C_3(:)-C_4(:)-C_2(:)-C_1(:))+W_6(:)-W_0(:))

        xdUM(:) = 0.125*(0.5*(W_4(:)+W_3(:)-C_4(:)-C_1(:)+W_5(:)+W_6(:)-E_5(:)-E_1(:))+ &
                    0.25*(E_3(:)+E_4(:)+E_7(:)+E_8(:)-C_5(:)-C_8(:)-C_9(:)-C_12(:))+W_0(:)-W_1(:))

        ydUM(:) = 0.125*(0.5*(W_2(:)+W_1(:)-C_2(:)-C_1(:)+W_6(:)+W_5(:)-E_3(:)-E_7(:))+ & 
                    0.25*(E_6(:)+E_5(:)+E_2(:)+E_1(:)-C_10(:)-C_9(:)-C_6(:)-C_5(:))+W_0(:)-W_3(:))

        zdUM(:) = 0.125*(0.5*(W_2(:)+W_1(:)-E_2(:)-E_1(:)+W_4(:)+W_3(:)-E_3(:)-E_4(:))+ & 
                    0.25*(C_3(:)+C_4(:)+C_2(:)+C_1(:)-C_7(:)-C_8(:)-C_6(:)-C_5(:))+W_0(:)-W_5(:))

        tmp_pos(1,:) = W_0(:)-(phi/4.)*((1+k)*xdUM(:)+(1-k)*xdUP(:))
        tmp_pos(2,:) = W_0(:)+(phi/4.)*((1-k)*xdUM(:)+(1+k)*xdUP(:))
        tmp_pos(3,:) = W_0(:)-(phi/4.)*((1+k)*ydUM(:)+(1-k)*ydUP(:))
        tmp_pos(4,:) = W_0(:)+(phi/4.)*((1-k)*ydUM(:)+(1+k)*ydUP(:))
        tmp_pos(5,:) = W_0(:)-(phi/4.)*((1+k)*zdUM(:)+(1-k)*zdUP(:))
        tmp_pos(6,:) = W_0(:)+(phi/4.)*((1-k)*zdUM(:)+(1+k)*zdUP(:))

        minpath(:) = min(W_0(:),W_1(:),W_2(:),W_3(:),W_4(:),W_5(:),W_6(:),C_1(:),C_2(:),C_3(:),C_4(:),C_5(:), &
                         C_6(:),C_7(:),C_8(:),C_9(:),C_10(:),C_11(:),C_12(:),E_1(:),E_2(:),E_3(:),E_4(:),E_5(:), &
                         E_6(:),E_7(:),E_8(:))
        maxpath(:) = max(W_0(:),W_1(:),W_2(:),W_3(:),W_4(:),W_5(:),W_6(:),C_1(:),C_2(:),C_3(:),C_4(:),C_5(:), &
                         C_6(:),C_7(:),C_8(:),C_9(:),C_10(:),C_11(:),C_12(:),E_1(:),E_2(:),E_3(:),E_4(:),E_5(:), &
                         E_6(:),E_7(:),E_8(:))

        dU2(:) = (maxpath(:)-minpath(:))**2. 
        limiter1(:) = 1.

        do i= 1, ele_nvar
           do face= 1, 6
              if ((tmp_pos(face,i)-W_0(i))>0.) then
                 y(i) = (maxpath(i)-W_0(i))/(tmp_pos(face,i)-W_0(i))
              end if
              if ((tmp_pos(face,i)-W_0(i))<0.) then
                 y(i) = (minpath(i)-W_0(i))/(tmp_pos(face,i)-W_0(i))
              end if
              if ((tmp_pos(face,i)-W_0(i))==0.) then
                 y(i) = 2.
              end if
              if (y(i)<1.5) then
                 poly(i) = -((4./27.)*(y(i)**3.))+y(i)
              else
                 poly(i) = 1.
              end if
              limiter1(i) = min(limiter1(i),poly(i))
           end do
        
           if (dU2(i)<trigger) then
              sig(i) = 1.
           end if
           if ((dU2(i)>=trigger).and.(dU2(i)<2.*trigger)) then
              sig(i) = 2.*(((dU2(i)-trigger)/trigger)**3.)-3.*(((dU2(i)-trigger)/trigger)**2.)+1.
           end if
           if (dU2(i)>=2.*trigger) then
              sig(i) = 0.
           end if
           limiter(i) = sig(i)+(1-sig(i))*limiter1(i)
        end do

        W_0_edge1(:) = W_0(:)-limiter(:)*(phi/4.)*((1+k)*xdUM(:)+(1-k)*xdUP(:))
        W_0_edge2(:) = W_0(:)+limiter(:)*(phi/4.)*((1-k)*xdUM(:)+(1+k)*xdUP(:))
        W_0_edge3(:) = W_0(:)-limiter(:)*(phi/4.)*((1+k)*ydUM(:)+(1-k)*ydUP(:))
        W_0_edge4(:) = W_0(:)+limiter(:)*(phi/4.)*((1-k)*ydUM(:)+(1+k)*ydUP(:))
        W_0_edge5(:) = W_0(:)-limiter(:)*(phi/4.)*((1+k)*zdUM(:)+(1-k)*zdUP(:))
        W_0_edge6(:) = W_0(:)+limiter(:)*(phi/4.)*((1-k)*zdUM(:)+(1+k)*zdUP(:))
  
        return

        end subroutine ele_reconstruct_edge

!-------------------------------------------------------------------------------!   
! Subroutine ele_reconstruct_edge1                                              !
!-------------------------------------------------------------------------------!

        subroutine ele_reconstruct_edge1(lb, W_0, W_1, W_2, W_3, W_4, W_5, W_6, &
        C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12, E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, &
        W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7, ord, dx)

        implicit none

        integer, intent(in) :: ord, lb
        real, intent(in) :: dx
        real, dimension(ele_nvar), intent(in) :: W_0, W_1, W_2, W_3, W_4, W_5, W_6
        real, dimension(ele_nvar), intent(in) :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12
        real, dimension(ele_nvar), intent(in) :: E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8
        real, dimension(ele_nvar), intent(out):: W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7
        real, dimension(ele_nvar) :: tmp_pos1, tmp_pos2, tmp_pos3, tmp_pos4, tmp_pos5, tmp_pos6, tmp_pos7, tmp_pos8
        real, dimension(ele_nvar) :: grad_x, grad_y, grad_z
        real, dimension(ele_nvar) ::  limiter, minpath, maxpath, mincell, maxcell
        real :: term1, term2, s, lim
        integer :: i
        
        W_0_edge7(:) = W_0(:)
        if (ord .eq. 1) then
           W_0_edge1(:) = W_0(:)
           W_0_edge2(:) = W_0(:)
           W_0_edge3(:) = W_0(:)
           W_0_edge4(:) = W_0(:)
           W_0_edge5(:) = W_0(:)
           W_0_edge6(:) = W_0(:)
           return
        end if

        s = 0.5

        grad_x(:) = 0.125*(0.5*(C_3(:)+C_2(:)-C_4(:)-C_1(:)+E_2(:)+E_6(:)-E_5(:)-E_1(:))+ &
                    0.25*(C_6(:)+C_7(:)+C_10(:)+C_11(:)-C_5(:)-C_8(:)-C_9(:)-C_12(:))+W_2(:)-W_1(:))

        grad_y(:) = 0.125*(0.5*(C_3(:)+C_4(:)-C_2(:)-C_1(:)+E_8(:)+E_4(:)-E_3(:)-E_7(:))+ & 
                    0.25*(C_11(:)+C_12(:)+C_7(:)+C_8(:)-C_10(:)-C_9(:)-C_6(:)-C_5(:))+W_4(:)-W_3(:))

        grad_z(:) = 0.125*(0.5*(E_6(:)+E_5(:)-E_2(:)-E_1(:)+E_8(:)+E_7(:)-E_3(:)-E_4(:))+ & 
                    0.25*(C_11(:)+C_12(:)+C_10(:)+C_9(:)-C_7(:)-C_8(:)-C_6(:)-C_5(:))+W_6(:)-W_5(:))


        minpath(:) = min(W_0(:),W_1(:),W_2(:),W_3(:),W_4(:),W_5(:),W_6(:),C_1(:),C_2(:),C_3(:),C_4(:),C_5(:), &
                         C_6(:),C_7(:),C_8(:),C_9(:),C_10(:),C_11(:),C_12(:),E_1(:),E_2(:),E_3(:),E_4(:),E_5(:), &
                         E_6(:),E_7(:),E_8(:))

        maxpath(:) = max(W_0(:),W_1(:),W_2(:),W_3(:),W_4(:),W_5(:),W_6(:),C_1(:),C_2(:),C_3(:),C_4(:),C_5(:), &
                         C_6(:),C_7(:),C_8(:),C_9(:),C_10(:),C_11(:),C_12(:),E_1(:),E_2(:),E_3(:),E_4(:),E_5(:), &
                         E_6(:),E_7(:),E_8(:))

        tmp_pos1(:) = W_0(:)-(s*grad_x(:))-(s*grad_y(:))-(s*grad_z(:))
        tmp_pos2(:) = W_0(:)+(s*grad_x(:))-(s*grad_y(:))-(s*grad_z(:))
        tmp_pos3(:) = W_0(:)+(s*grad_x(:))+(s*grad_y(:))-(s*grad_z(:))
        tmp_pos4(:) = W_0(:)-(s*grad_x(:))+(s*grad_y(:))-(s*grad_z(:))
        tmp_pos5(:) = W_0(:)-(s*grad_x(:))-(s*grad_y(:))+(s*grad_z(:))
        tmp_pos6(:) = W_0(:)+(s*grad_x(:))-(s*grad_y(:))+(s*grad_z(:))
        tmp_pos7(:) = W_0(:)+(s*grad_x(:))+(s*grad_y(:))+(s*grad_z(:))
        tmp_pos8(:) = W_0(:)-(s*grad_x(:))+(s*grad_y(:))+(s*grad_z(:))

        mincell(:) = min(tmp_pos1(:),tmp_pos2(:),tmp_pos3(:),tmp_pos4(:),tmp_pos5(:),tmp_pos6(:),tmp_pos7(:),tmp_pos8(:))
        maxcell(:) = max(tmp_pos1(:),tmp_pos2(:),tmp_pos3(:),tmp_pos4(:),tmp_pos5(:),tmp_pos6(:),tmp_pos7(:),tmp_pos8(:))
        
        do i=1,ele_nvar
           if (W_0(i) .ne. maxcell(i)) then 
              term1 = abs(W_0(i)-maxpath(i))/abs(W_0(i)-maxcell(i))
           else
              term1 = 1.
           end if
           if (W_0(i) .ne. mincell(i)) then 
              term2 = abs(W_0(i)-minpath(i))/abs(W_0(i)-mincell(i))
           else
              term2 = 1.
           end if
           
           limiter(i) = min(1.,term1, term2)
        end do    

        W_0_edge1(:) = W_0(:)-(0.5*limiter(:)*grad_x(:))
        W_0_edge2(:) = W_0(:)+(0.5*limiter(:)*grad_x(:))
        W_0_edge3(:) = W_0(:)-(0.5*limiter(:)*grad_y(:))
        W_0_edge4(:) = W_0(:)+(0.5*limiter(:)*grad_y(:))
        W_0_edge5(:) = W_0(:)-(0.5*limiter(:)*grad_z(:))
        W_0_edge6(:) = W_0(:)+(0.5*limiter(:)*grad_z(:))

        return

        end subroutine ele_reconstruct_edge1

!-------------------------------------------------------------------------------!   
! Subroutine ele_reconstruct_edge2                                              !
!-------------------------------------------------------------------------------!

        subroutine ele_reconstruct_edge2(lb, ic, jc, kc, W_0, W_1, W_2, W_3, W_4, W_5, W_6, &
        C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12, E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, &
        W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7, ord)


        implicit none

        integer, intent(in) :: ord, lb, ic, jc, kc
        real, dimension(ele_nvar), intent(in) :: W_0, W_1, W_2, W_3, W_4, W_5, W_6
        real, dimension(ele_nvar), intent(in) :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12
        real, dimension(ele_nvar), intent(in) :: E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8
        real, dimension(ele_nvar), intent(out):: W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7
        real, dimension(26) :: dx, dy, dz
        real, dimension(ele_nvar,-1:1,-1:1,-1:1) :: Wdat
        real, dimension(ele_nvar,26) :: dWdat
        real, dimension(ele_nvar) ::  sum_Wx, sum_Wy, sum_Wz, limiter, minpath, maxpath, mincell, maxcell
        real :: sum_x2, sum_y2, sum_z2, sum_xy, sum_xz, sum_yz, Den
        real, dimension(ele_nvar) :: tmp_pos1, tmp_pos2, tmp_pos3, tmp_pos4, tmp_pos5, tmp_pos6, tmp_pos7, tmp_pos8
        real, dimension(ele_nvar) :: grad_x, grad_y, grad_z
        real, dimension(7) :: fx, fy, fz
        real :: term1, term2, s, dsub, xc, yc, zc, xw, yw, zw, xref, yref, zref, xmin, xmax, ymin, ymax, zmin, zmax
        integer :: i, j, k, l, index, face

        if (ord .eq. 1) then
           W_0_edge1(:) = W_0(:)
           W_0_edge2(:) = W_0(:)
           W_0_edge3(:) = W_0(:)
           W_0_edge4(:) = W_0(:)
           W_0_edge5(:) = W_0(:)
           W_0_edge6(:) = W_0(:)
           W_0_edge7(:) = W_0(:)
           return
        end if

        dsub = real(bsize(1,lb)/nxb)
        s = 0.5*dsub
        xref = coord(1,lb)+(((2*(ic-nguard))-nxb-1)*bsize(1,lb)/(2*nxb))
        yref = coord(2,lb)+(((2*(jc-nguard))-nyb-1)*bsize(2,lb)/(2*nyb))
        zref = coord(3,lb)+(((2*(kc-nguard))-nzb-1)*bsize(3,lb)/(2*nzb))
        xmin = xref-s 
        xmax = xref+s 
        ymin = yref-s 
        ymax = yref+s 
        zmin = zref-s 
        zmax = zref+s 

        if (loc_cell_flags(ic,jc,kc).ne.1) then
           xc = xref
           yc = yref
           zc = zref
        else
           xc = surf_cell_data(2,ic,jc,kc,bflags(1,lb))
           yc = surf_cell_data(3,ic,jc,kc,bflags(1,lb))
           zc = surf_cell_data(4,ic,jc,kc,bflags(1,lb))
        end if
        
        Wdat(:,-1,-1,-1) = C_5(:)
        Wdat(:, 0,-1,-1) = E_3(:)
        Wdat(:, 1,-1,-1) = C_6(:)
        Wdat(:,-1, 0,-1) = E_1(:)
        Wdat(:, 0, 0,-1) = W_5(:)
        Wdat(:, 1, 0,-1) = E_2(:)
        Wdat(:,-1, 1,-1) = C_8(:)
        Wdat(:, 0, 1,-1) = E_4(:)
        Wdat(:, 1, 1,-1) = C_7(:)
        Wdat(:,-1,-1, 0) = C_1(:)
        Wdat(:, 0,-1, 0) = W_3(:)
        Wdat(:, 1,-1, 0) = C_2(:)
        Wdat(:,-1, 0, 0) = W_1(:)
        Wdat(:, 0, 0, 0) = W_0(:)
        Wdat(:, 1, 0, 0) = W_2(:)
        Wdat(:,-1, 1, 0) = C_4(:)
        Wdat(:, 0, 1, 0) = W_4(:)
        Wdat(:, 1, 1, 0) = C_3(:)
        Wdat(:,-1,-1, 1) = C_9(:)
        Wdat(:, 0,-1, 1) = E_7(:)
        Wdat(:, 1,-1, 1) = C_10(:)
        Wdat(:,-1, 0, 1) = E_5(:)
        Wdat(:, 0, 0, 1) = W_6(:)
        Wdat(:, 1, 0, 1) = E_6(:)
        Wdat(:,-1, 1, 1) = C_12(:)
        Wdat(:, 0, 1, 1) = E_8(:)
        Wdat(:, 1, 1, 1) = C_11(:)

        index = 0
        minpath(:) = Wdat(:,0,0,0)
        maxpath(:) = Wdat(:,0,0,0)
        do k = -1,1
           do j = -1,1
              do i= -1,1
                 if (((i.ne.0).or.(j.ne.0).or.(k.ne.0)).and.(loc_cell_flags(ic+i,jc+j,kc+k).ge.0)) then
                    index = index+1
                    if (loc_cell_flags(ic+i,jc+j,kc+k).eq.1) then
                       xw = surf_cell_data(2,ic+i,jc+j,kc+k,bflags(1,lb))
                       yw = surf_cell_data(3,ic+i,jc+j,kc+k,bflags(1,lb))
                       zw = surf_cell_data(4,ic+i,jc+j,kc+k,bflags(1,lb))
                    else
                       xw = xref +(real(i)*dsub)
                       yw = yref +(real(j)*dsub)
                       zw = zref +(real(k)*dsub)
                    end if
                    dWdat(:,index) = Wdat(:,i,j,k)-Wdat(:,0,0,0)
                    dx(index) = xw-xc
                    dy(index) = yw-yc
                    dz(index) = zw-zc
                    minpath(:) =min(minpath(:),Wdat(:,i,j,k))
                    maxpath(:) =max(maxpath(:),Wdat(:,i,j,k))
                 end if
              end do
           end do
        end do

        sum_Wx(:) = 0.
        sum_Wy(:) = 0.
        sum_Wz(:) = 0.
        sum_x2 = 0.
        sum_y2 = 0.
        sum_z2 = 0.
        sum_xy = 0.
        sum_xz = 0.
        sum_yz = 0.

        do l=1,index
           sum_Wx(:) = sum_Wx(:)+(dWdat(:,l)*dx(l))
           sum_Wy(:) = sum_Wy(:)+(dWdat(:,l)*dy(l))
           sum_Wz(:) = sum_Wz(:)+(dWdat(:,l)*dz(l))
           sum_x2 = sum_x2+(dx(l)**2.)
           sum_y2 = sum_y2+(dy(l)**2.)
           sum_z2 = sum_z2+(dz(l)**2.)
           sum_xy = sum_xy+(dx(l)*dy(l))
           sum_xz = sum_xz+(dx(l)*dz(l))
           sum_yz = sum_yz+(dy(l)*dz(l))
        end do

        Den=(sum_x2*sum_y2*sum_z2)-sum_x2*(sum_yz**2.)-sum_y2*(sum_xz**2.)-sum_z2*(sum_xy**2.)+(2.*sum_xy*sum_xz*sum_yz)

        grad_x(:)=(((sum_y2*sum_z2-(sum_yz**2.))*sum_Wx(:))+((sum_xz*sum_yz-sum_xy*sum_z2)*sum_Wy(:))+ &
                    ((sum_xy*sum_yz-sum_xz*sum_y2)*sum_Wz(:))/Den)
        grad_y(:)=(((sum_xz*sum_yz-sum_xy*sum_z2)*sum_Wx(:))+((sum_x2*sum_z2-(sum_xz**2.))*sum_Wy(:))+ &
                   ((sum_xy*sum_xz-sum_yz*sum_x2)*sum_Wz(:))/Den)
        grad_z(:)=(((sum_xy*sum_yz-sum_xz*sum_y2)*sum_Wx(:))+((sum_xy*sum_xz-sum_yz*sum_x2)*sum_Wy(:))+ &
                   ((sum_x2*sum_y2-(sum_xy**2.))*sum_Wz(:))/Den)

        mincell(:) = minpath(:)
        maxcell(:) = maxpath(:)
        do face=1,7
           fx(face) = surf_face_data(8,face,ic,jc,kc,bflags(1,lb))-xc
           fy(face) = surf_face_data(9,face,ic,jc,kc,bflags(1,lb))-yc
           fz(face) = surf_face_data(10,face,ic,jc,kc,bflags(1,lb))-zc
           if (surf_face_data(3,face,ic,jc,kc,bflags(1,lb)).ne.0) then
              tmp_pos1(:) = W_0(:)+(fx(face)*grad_x(:))+(fy(face)*grad_y(:))+(fz(face)*grad_z(:))
              mincell(:) = min(mincell(:),tmp_pos1(:))
              maxcell(:) = max(maxcell(:),tmp_pos1(:))
           end if
        end do
        
        do i=1,ele_nvar
           if (W_0(i) .ne. maxcell(i)) then 
              term1 = abs(W_0(i)-maxpath(i))/abs(W_0(i)-maxcell(i))
           else
              term1 = 1.
           end if
           if (W_0(i) .ne. mincell(i)) then 
              term2 = abs(W_0(i)-minpath(i))/abs(W_0(i)-mincell(i))
           else
              term2 = 1.
           end if
           
           limiter(i) = min(1.,term1, term2)
        end do   
        
        W_0_edge1(:) = W_0(:)+(fx(1)*limiter(:)*grad_x(:))+(fy(1)*limiter(:)*grad_y(:))+(fz(1)*limiter(:)*grad_z(:))
        W_0_edge2(:) = W_0(:)+(fx(2)*limiter(:)*grad_x(:))+(fy(2)*limiter(:)*grad_y(:))+(fz(2)*limiter(:)*grad_z(:))
        W_0_edge3(:) = W_0(:)+(fx(3)*limiter(:)*grad_x(:))+(fy(3)*limiter(:)*grad_y(:))+(fz(3)*limiter(:)*grad_z(:))
        W_0_edge4(:) = W_0(:)+(fx(4)*limiter(:)*grad_x(:))+(fy(4)*limiter(:)*grad_y(:))+(fz(4)*limiter(:)*grad_z(:))
        W_0_edge5(:) = W_0(:)+(fx(5)*limiter(:)*grad_x(:))+(fy(5)*limiter(:)*grad_y(:))+(fz(5)*limiter(:)*grad_z(:))
        W_0_edge6(:) = W_0(:)+(fx(6)*limiter(:)*grad_x(:))+(fy(6)*limiter(:)*grad_y(:))+(fz(6)*limiter(:)*grad_z(:))
        W_0_edge7(:) = W_0(:)+(fx(7)*limiter(:)*grad_x(:))+(fy(7)*limiter(:)*grad_y(:))+(fz(7)*limiter(:)*grad_z(:))

        return

        end subroutine ele_reconstruct_edge2

end module ele_module
         
