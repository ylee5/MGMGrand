
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
!                                FIELD_MODULE.F90                               !
!-------------------------------------------------------------------------------!
! Header 
        module field_module

        use declarationmodule
        use simulationparameters
        use misc_module
        use ieee_arithmetic

        implicit none

        contains 

!-------------------------------------------------------------------------------!   
! Subroutine field_face_normal                                                  !
!-------------------------------------------------------------------------------! 

        subroutine field_face_normal(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k

! Interpolate centred values to face values -----------------------------------!
        call  field_interpol_to_faces_normal(lb, istep, G0_B_prim(:,:,:,:,lb), B_face_xm(:,:,:,:,lb), B_face_xp(:,:,:,:,lb), &
              B_face_ym(:,:,:,:,lb), B_face_yp(:,:,:,:,lb), B_face_zm(:,:,:,:,lb), B_face_zp(:,:,:,:,lb), order)

! Compute the average face values ---------------------------------------------!
        B_face_xm_av(:,:,:,:,lb)=0.
        B_face_xp_av(:,:,:,:,lb)=0.
        B_face_ym_av(:,:,:,:,lb)=0.
        B_face_yp_av(:,:,:,:,lb)=0.
        B_face_zm_av(:,:,:,:,lb)=0.

        do i=nguard+1,nxb+nguard
           do j=nguard+1,nyb+nguard
              do k=nguard+1,nzb+nguard
                 B_face_xm_av(:,i,j,k,lb)=0.5*(B_face_xm(:,i,j,k,lb)+B_face_xp(:,i-1,j,k,lb))
                 B_face_xp_av(:,i,j,k,lb)=0.5*(B_face_xp(:,i,j,k,lb)+B_face_xm(:,i+1,j,k,lb))
                 B_face_ym_av(:,i,j,k,lb)=0.5*(B_face_ym(:,i,j,k,lb)+B_face_yp(:,i,j-1,k,lb))
                 B_face_yp_av(:,i,j,k,lb)=0.5*(B_face_yp(:,i,j,k,lb)+B_face_ym(:,i,j+1,k,lb))
                 B_face_zm_av(:,i,j,k,lb)=0.5*(B_face_zm(:,i,j,k,lb)+B_face_zp(:,i,j,k-1,lb))
                 B_face_zp_av(:,i,j,k,lb)=0.5*(B_face_zp(:,i,j,k,lb)+B_face_zm(:,i,j,k+1,lb))
              end do
           end do
        end do

        end subroutine field_face_normal

!-------------------------------------------------------------------------------!   
! Subroutine field_face_surface                                                 !
!-------------------------------------------------------------------------------! 

        subroutine field_face_surface(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k

! Interpolate centred values to face values -----------------------------------!
        call  field_interpol_to_faces_surface(lb, istep, G0_B_prim(:,:,:,:,lb), B_face_xm(:,:,:,:,lb), B_face_xp(:,:,:,:,lb), &
              B_face_ym(:,:,:,:,lb), B_face_yp(:,:,:,:,lb), B_face_zm(:,:,:,:,lb), B_face_zp(:,:,:,:,lb), B_face_ei(:,:,:,:,lb),&
              B_face_eo(:,:,:,:,lb),order)

! Compute the average face values ---------------------------------------------!
        B_face_xm_av(:,:,:,:,lb)=0.
        B_face_xp_av(:,:,:,:,lb)=0.
        B_face_ym_av(:,:,:,:,lb)=0.
        B_face_yp_av(:,:,:,:,lb)=0.
        B_face_zm_av(:,:,:,:,lb)=0.
        B_face_zp_av(:,:,:,:,lb)=0.
        B_face_ex_av(:,:,:,:,lb)=0.

        do i=nguard+1,nxb+nguard
           do j=nguard+1,nyb+nguard
              do k=nguard+1,nzb+nguard
                 B_face_xm_av(:,i,j,k,lb)=0.5*(B_face_xm(:,i,j,k,lb)+B_face_xp(:,i-1,j,k,lb))
                 B_face_xp_av(:,i,j,k,lb)=0.5*(B_face_xp(:,i,j,k,lb)+B_face_xm(:,i+1,j,k,lb))
                 B_face_ym_av(:,i,j,k,lb)=0.5*(B_face_ym(:,i,j,k,lb)+B_face_yp(:,i,j-1,k,lb))
                 B_face_yp_av(:,i,j,k,lb)=0.5*(B_face_yp(:,i,j,k,lb)+B_face_ym(:,i,j+1,k,lb))
                 B_face_zm_av(:,i,j,k,lb)=0.5*(B_face_zm(:,i,j,k,lb)+B_face_zp(:,i,j,k-1,lb))
                 B_face_zp_av(:,i,j,k,lb)=0.5*(B_face_zp(:,i,j,k,lb)+B_face_zm(:,i,j,k+1,lb))
                 B_face_ex_av(:,i,j,k,lb)=0.5*(B_face_ei(:,i,j,k,lb)+B_face_eo(:,i,j,k,lb))
              end do
           end do
        end do

        end subroutine field_face_surface

!-------------------------------------------------------------------------------!   
! Subroutine current_face_normal                                                !
!-------------------------------------------------------------------------------! 

        subroutine current_face_normal(mype, lb)

        implicit none 

        integer, intent(in) ::lb, mype
        integer :: i, j, k, ii, jj, kk
        real, dimension(3) :: B1, B2, B3, B4, B5, B6
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: B
        real :: dx, ele_nd
        real, dimension(3) :: rot_B, Jd

        dx = (bsize(1,lb)*length_unit/(l_sc*nxb))
 
! Compute the current field ----------------------------------------------------!

        if (hall_mhd_flag==.true.) then
           
           B(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) = (G0_B_prim(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard,lb)+B0_prim(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard,lb))*B_sc
           
           do k=1,nzb                                                        
              kk = k+nguard
              do j=1,nyb
                 jj = j+nguard
                 do i=1,nxb
                    ii= i+nguard
                    
                    J_face_xp(1,ii,jj,kk,lb) = (B(3,ii,jj+1,kk)+B(3,ii+1,jj+1,kk)-B(3,ii,jj-1,kk)-B(3,ii+1,jj-1,kk)-B(2,ii,jj,kk+1)-B(2,ii+1,jj,kk+1)+B(2,ii,jj,kk-1)+B(2,ii+1,jj,kk-1))/(4.*dx*l_sc*mu_0)
                    J_face_xp(2,ii,jj,kk,lb) = (B(1,ii,jj,kk+1)+B(1,ii+1,jj,kk+1)-B(1,ii,jj,kk-1)-B(1,ii+1,jj,kk-1)-4.*B(3,ii+1,jj,kk)+4.*B(3,ii,jj,kk))/(4.*dx*l_sc*mu_0)
                    J_face_xp(3,ii,jj,kk,lb) = (4.*B(2,ii+1,jj,kk)-4.*B(2,ii,jj,kk)-B(1,ii,jj+1,kk)-B(1,ii+1,jj+1,kk)+B(1,ii,jj-1,kk)+B(1,ii+1,jj-1,kk))/(4.*dx*l_sc*mu_0)
                    J_face_xm(1,ii,jj,kk,lb) = (B(3,ii-1,jj+1,kk)+B(3,ii,jj+1,kk)-B(3,ii-1,jj-1,kk)-B(3,ii,jj-1,kk)-B(2,ii-1,jj,kk+1)-B(2,ii,jj,kk+1)+B(2,ii-1,jj,kk-1)+B(2,ii,jj,kk-1))/(4.*dx*l_sc*mu_0)
                    J_face_xm(2,ii,jj,kk,lb) = (B(1,ii-1,jj,kk+1)+B(1,ii,jj,kk+1)-B(1,ii-1,jj,kk-1)-B(1,ii,jj,kk-1)-4.*B(3,ii,jj,kk)+4.*B(3,ii-1,jj,kk))/(4.*dx*l_sc*mu_0)
                    J_face_xm(3,ii,jj,kk,lb) = (4.*B(2,ii,jj,kk)-4.*B(2,ii-1,jj,kk)-B(1,ii-1,jj+1,kk)-B(1,ii,jj+1,kk)+B(1,ii-1,jj-1,kk)+B(1,ii,jj-1,kk))/(4.*dx*l_sc*mu_0)  
                    J_face_yp(1,ii,jj,kk,lb) = (B(2,ii,jj,kk+1)+B(2,ii,jj+1,kk+1)-B(2,ii,jj,kk-1)-B(2,ii,jj+1,kk-1)-4.*B(3,ii,jj+1,kk)+4.*B(3,ii,jj,kk))/(4.*dx*l_sc*mu_0)  
                    J_face_yp(2,ii,jj,kk,lb) = (B(3,ii+1,jj,kk)+B(3,ii+1,jj+1,kk)-B(3,ii-1,jj,kk)-B(3,ii-1,jj+1,kk)-B(1,ii,jj,kk+1)-B(1,ii,jj+1,kk+1)+B(1,ii,jj,kk-1)+B(1,ii,jj+1,kk-1))/(4.*dx*l_sc*mu_0)
                    J_face_yp(3,ii,jj,kk,lb) = (4.*B(1,ii,jj+1,kk)-4.*B(1,ii,jj,kk)-B(2,ii+1,jj,kk)-B(2,ii+1,jj+1,kk)+B(2,ii-1,jj,kk)+B(2,ii-1,jj+1,kk))/(4.*dx*l_sc*mu_0)
                    J_face_ym(1,ii,jj,kk,lb) = (B(2,ii,jj-1,kk+1)+B(2,ii,jj,kk+1)-B(2,ii,jj-1,kk-1)-B(2,ii,jj,kk-1)-4.*B(3,ii,jj,kk)+4.*B(3,ii,jj-1,kk))/(4.*dx*l_sc*mu_0)
                    J_face_ym(2,ii,jj,kk,lb) = (B(3,ii+1,jj-1,kk)+B(3,ii+1,jj,kk)-B(3,ii-1,jj-1,kk)-B(3,ii-1,jj,kk)-B(1,ii,jj-1,kk+1)-B(1,ii,jj,kk+1)+B(1,ii,jj-1,kk-1)+B(1,ii,jj,kk-1))/(4.*dx*l_sc*mu_0)
                    J_face_ym(3,ii,jj,kk,lb) = (4.*B(1,ii,jj,kk)-4.*B(1,ii,jj-1,kk)-B(2,ii+1,jj-1,kk)-B(2,ii+1,jj,kk)+B(2,ii-1,jj-1,kk)+B(2,ii-1,jj,kk))/(4.*dx*l_sc*mu_0)
                    J_face_zp(1,ii,jj,kk,lb) = (4.*B(2,ii,jj,kk+1)-4.*B(2,ii,jj,kk)-B(3,ii,jj+1,kk)-B(3,ii,jj+1,kk+1)+B(3,ii,jj-1,kk)+B(3,ii,jj-1,kk+1))/(4.*dx*l_sc*mu_0)
                    J_face_zp(2,ii,jj,kk,lb) = (B(3,ii+1,jj,kk)+B(3,ii+1,jj,kk+1)-B(3,ii-1,jj,kk)-B(3,ii-1,jj,kk+1)-4.*B(1,ii,jj,kk+1)+4.*B(1,ii,jj,kk))/(4.*dx*l_sc*mu_0)
                    J_face_zp(3,ii,jj,kk,lb) = (B(1,ii,jj+1,kk)+B(1,ii,jj+1,kk+1)-B(1,ii,jj-1,kk)-B(1,ii,jj-1,kk+1)-B(2,ii+1,jj,kk)-B(2,ii+1,jj,kk+1)+B(2,ii-1,jj,kk)+B(2,ii-1,jj,kk+1))/(4.*dx*l_sc*mu_0)
                    J_face_zm(1,ii,jj,kk,lb) = (4.*B(2,ii,jj,kk)-4.*B(2,ii,jj,kk-1)-B(3,ii,jj+1,kk-1)-B(3,ii,jj+1,kk)+B(3,ii,jj-1,kk-1)+B(3,ii,jj-1,kk))/(4.*dx*l_sc*mu_0)
                    J_face_zm(2,ii,jj,kk,lb) = (B(3,ii+1,jj,kk-1)+B(3,ii+1,jj,kk)-B(3,ii-1,jj,kk-1)-B(3,ii-1,jj,kk)-4.*B(1,ii,jj,kk)+4.*B(1,ii,jj,kk-1))/(4.*dx*l_sc*mu_0)
                    J_face_zm(3,ii,jj,kk,lb) = (B(1,ii,jj+1,kk-1)+B(1,ii,jj+1,kk)-B(1,ii,jj-1,kk-1)-B(1,ii,jj-1,kk)-B(2,ii+1,jj,kk-1)-B(2,ii+1,jj,kk)+B(2,ii-1,jj,kk-1)+B(2,ii-1,jj,kk))/(4.*dx*l_sc*mu_0)  

                    B1(:) = (B_face_xm(:,ii,jj,kk,lb)+B0_face_xm(:,ii,jj,kk,lb))*B_sc
                    B2(:) = (B_face_xp(:,ii,jj,kk,lb)+B0_face_xp(:,ii,jj,kk,lb))*B_sc
                    B3(:) = (B_face_ym(:,ii,jj,kk,lb)+B0_face_ym(:,ii,jj,kk,lb))*B_sc
                    B4(:) = (B_face_yp(:,ii,jj,kk,lb)+B0_face_yp(:,ii,jj,kk,lb))*B_sc
                    B5(:) = (B_face_zm(:,ii,jj,kk,lb)+B0_face_zm(:,ii,jj,kk,lb))*B_sc
                    B6(:) = (B_face_zp(:,ii,jj,kk,lb)+B0_face_zp(:,ii,jj,kk,lb))*B_sc

                    ele_nd =(G0_ion_prim(1,ii,jj,kk,lb)*rho_sc)/(M_p*amu)
                    rot_B(:) =  (/B4(3)-B3(3)+B5(2)-B6(2), B1(3)-B2(3)+B6(1)-B5(1), B2(2)-B1(2)+B3(1)-B4(1)/) 
                    Jd(:) = Rot_B(:)/(dx*l_sc*mu_0)
                    G0_ele_prim(2:4,ii,jj,kk,lb) = G0_ion_prim(2:4,ii,jj,kk,lb)-(Jd(:)/(ele_nd*q*a_sc))
                    J_prim(:,ii,jj,kk,lb) = Jd(:)/J_sc
                    
                 end do
              end do
           end do
           
        else

           J_face_xp(:,:,:,:,lb) = 0.
           J_face_xm(:,:,:,:,lb) = 0.
           J_face_yp(:,:,:,:,lb) = 0.
           J_face_ym(:,:,:,:,lb) = 0.
           J_face_zp(:,:,:,:,lb) = 0.
           J_face_zm(:,:,:,:,lb) = 0.
           G0_ele_prim(2:4,:,:,:,lb) = G0_ion_prim(2:4,:,:,:,lb)
           J_prim(:,:,:,:,lb) = 0.

        end if

        G0_ele_prim(1,:,:,:,lb) = (G0_ion_prim(1,:,:,:,lb)*rho_sc)/(M_p*amu*nd_sc)
        G0_ele_cons(1,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)
        G0_ele_cons(2,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)*G0_ele_prim(2,:,:,:,lb)
        G0_ele_cons(3,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)*G0_ele_prim(3,:,:,:,lb)
        G0_ele_cons(4,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)*G0_ele_prim(4,:,:,:,lb)
        G0_ele_cons(5,:,:,:,lb) = G0_ele_prim(5,:,:,:,lb) 
        J_save(1:3,1:nxb,1:nyb,1:nzb,lb) = J_prim(1:3,nguard+1:nxb+nguard,nguard+1:nyb+nguard,nguard+1:nzb+nguard,lb)

        return

        end subroutine current_face_normal

!-------------------------------------------------------------------------------!   
! Subroutine current_face_surface                                               !
!-------------------------------------------------------------------------------! 

        subroutine current_face_surface(mype, lb)

        implicit none 

        integer, intent(in) ::lb, mype
        integer :: i, j, k, ii, jj, kk, face
        real, dimension(3) :: B1, B2, B3, B4, B5, B6, B7
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: B
        real :: nx, ny, nz
        real :: dx, dx_scale, volume, surface, ele_nd
        real, dimension(3) :: rot_B, Jd

        dx = (bsize(1,lb)*length_unit/(l_sc*nxb))
        dx_scale = length_unit/l_sc

! Compute the current field ----------------------------------------------------!

        if (hall_mhd_flag==.true.) then

           B(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) = (G0_B_prim(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard,lb)+B0_prim(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard,lb))*B_sc

           do k=1,nzb   
              kk = k+nguard                            
              do j=1,nyb
                 jj = j+nguard
                 do i=1,nxb
                    ii= i+nguard

! Case of a non surface cell ----------------------------------------------------!
                    if (loc_cell_flags(ii,jj,kk).eq.0) then
            
                       J_face_xp(1,ii,jj,kk,lb) = (B(3,ii,jj+1,kk)+B(3,ii+1,jj+1,kk)-B(3,ii,jj-1,kk)-B(3,ii+1,jj-1,kk)-B(2,ii,jj,kk+1)-B(2,ii+1,jj,kk+1)+B(2,ii,jj,kk-1)+B(2,ii+1,jj,kk-1))/(4.*dx*l_sc*mu_0)
                       J_face_xp(2,ii,jj,kk,lb) = (B(1,ii,jj,kk+1)+B(1,ii+1,jj,kk+1)-B(1,ii,jj,kk-1)-B(1,ii+1,jj,kk-1)-4.*B(3,ii+1,jj,kk)+4.*B(3,ii,jj,kk))/(4.*dx*l_sc*mu_0)
                       J_face_xp(3,ii,jj,kk,lb) = (4.*B(2,ii+1,jj,kk)-4.*B(2,ii,jj,kk)-B(1,ii,jj+1,kk)-B(1,ii+1,jj+1,kk)+B(1,ii,jj-1,kk)+B(1,ii+1,jj-1,kk))/(4.*dx*l_sc*mu_0)
                       J_face_xm(1,ii,jj,kk,lb) = (B(3,ii-1,jj+1,kk)+B(3,ii,jj+1,kk)-B(3,ii-1,jj-1,kk)-B(3,ii,jj-1,kk)-B(2,ii-1,jj,kk+1)-B(2,ii,jj,kk+1)+B(2,ii-1,jj,kk-1)+B(2,ii,jj,kk-1))/(4.*dx*l_sc*mu_0)
                       J_face_xm(2,ii,jj,kk,lb) = (B(1,ii-1,jj,kk+1)+B(1,ii,jj,kk+1)-B(1,ii-1,jj,kk-1)-B(1,ii,jj,kk-1)-4.*B(3,ii,jj,kk)+4.*B(3,ii-1,jj,kk))/(4.*dx*l_sc*mu_0)
                       J_face_xm(3,ii,jj,kk,lb) = (4.*B(2,ii,jj,kk)-4.*B(2,ii-1,jj,kk)-B(1,ii-1,jj+1,kk)-B(1,ii,jj+1,kk)+B(1,ii-1,jj-1,kk)+B(1,ii,jj-1,kk))/(4.*dx*l_sc*mu_0)  
                       J_face_yp(1,ii,jj,kk,lb) = (B(2,ii,jj,kk+1)+B(2,ii,jj+1,kk+1)-B(2,ii,jj,kk-1)-B(2,ii,jj+1,kk-1)-4.*B(3,ii,jj+1,kk)+4.*B(3,ii,jj,kk))/(4.*dx*l_sc*mu_0)  
                       J_face_yp(2,ii,jj,kk,lb) = (B(3,ii+1,jj,kk)+B(3,ii+1,jj+1,kk)-B(3,ii-1,jj,kk)-B(3,ii-1,jj+1,kk)-B(1,ii,jj,kk+1)-B(1,ii,jj+1,kk+1)+B(1,ii,jj,kk-1)+B(1,ii,jj+1,kk-1))/(4.*dx*l_sc*mu_0)
                       J_face_yp(3,ii,jj,kk,lb) = (4.*B(1,ii,jj+1,kk)-4.*B(1,ii,jj,kk)-B(2,ii+1,jj,kk)-B(2,ii+1,jj+1,kk)+B(2,ii-1,jj,kk)+B(2,ii-1,jj+1,kk))/(4.*dx*l_sc*mu_0)
                       J_face_ym(1,ii,jj,kk,lb) = (B(2,ii,jj-1,kk+1)+B(2,ii,jj,kk+1)-B(2,ii,jj-1,kk-1)-B(2,ii,jj,kk-1)-4.*B(3,ii,jj,kk)+4.*B(3,ii,jj-1,kk))/(4.*dx*l_sc*mu_0)
                       J_face_ym(2,ii,jj,kk,lb) = (B(3,ii+1,jj-1,kk)+B(3,ii+1,jj,kk)-B(3,ii-1,jj-1,kk)-B(3,ii-1,jj,kk)-B(1,ii,jj-1,kk+1)-B(1,ii,jj,kk+1)+B(1,ii,jj-1,kk-1)+B(1,ii,jj,kk-1))/(4.*dx*l_sc*mu_0)
                       J_face_ym(3,ii,jj,kk,lb) = (4.*B(1,ii,jj,kk)-4.*B(1,ii,jj-1,kk)-B(2,ii+1,jj-1,kk)-B(2,ii+1,jj,kk)+B(2,ii-1,jj-1,kk)+B(2,ii-1,jj,kk))/(4.*dx*l_sc*mu_0)
                       J_face_zp(1,ii,jj,kk,lb) = (4.*B(2,ii,jj,kk+1)-4.*B(2,ii,jj,kk)-B(3,ii,jj+1,kk)-B(3,ii,jj+1,kk+1)+B(3,ii,jj-1,kk)+B(3,ii,jj-1,kk+1))/(4.*dx*l_sc*mu_0)
                       J_face_zp(2,ii,jj,kk,lb) = (B(3,ii+1,jj,kk)+B(3,ii+1,jj,kk+1)-B(3,ii-1,jj,kk)-B(3,ii-1,jj,kk+1)-4.*B(1,ii,jj,kk+1)+4.*B(1,ii,jj,kk))/(4.*dx*l_sc*mu_0)
                       J_face_zp(3,ii,jj,kk,lb) = (B(1,ii,jj+1,kk)+B(1,ii,jj+1,kk+1)-B(1,ii,jj-1,kk)-B(1,ii,jj-1,kk+1)-B(2,ii+1,jj,kk)-B(2,ii+1,jj,kk+1)+B(2,ii-1,jj,kk)+B(2,ii-1,jj,kk+1))/(4.*dx*l_sc*mu_0)
                       J_face_zm(1,ii,jj,kk,lb) = (4.*B(2,ii,jj,kk)-4.*B(2,ii,jj,kk-1)-B(3,ii,jj+1,kk-1)-B(3,ii,jj+1,kk)+B(3,ii,jj-1,kk-1)+B(3,ii,jj-1,kk))/(4.*dx*l_sc*mu_0)
                       J_face_zm(2,ii,jj,kk,lb) = (B(3,ii+1,jj,kk-1)+B(3,ii+1,jj,kk)-B(3,ii-1,jj,kk-1)-B(3,ii-1,jj,kk)-4.*B(1,ii,jj,kk)+4.*B(1,ii,jj,kk-1))/(4.*dx*l_sc*mu_0)
                       J_face_zm(3,ii,jj,kk,lb) = (B(1,ii,jj+1,kk-1)+B(1,ii,jj+1,kk)-B(1,ii,jj-1,kk-1)-B(1,ii,jj-1,kk)-B(2,ii+1,jj,kk-1)-B(2,ii+1,jj,kk)+B(2,ii-1,jj,kk-1)+B(2,ii-1,jj,kk))/(4.*dx*l_sc*mu_0)

                       B1(:) = (B_face_xm(:,ii,jj,kk,lb)+B0_face_xm(:,ii,jj,kk,lb))*B_sc
                       B2(:) = (B_face_xp(:,ii,jj,kk,lb)+B0_face_xp(:,ii,jj,kk,lb))*B_sc
                       B3(:) = (B_face_ym(:,ii,jj,kk,lb)+B0_face_ym(:,ii,jj,kk,lb))*B_sc
                       B4(:) = (B_face_yp(:,ii,jj,kk,lb)+B0_face_yp(:,ii,jj,kk,lb))*B_sc
                       B5(:) = (B_face_zm(:,ii,jj,kk,lb)+B0_face_zm(:,ii,jj,kk,lb))*B_sc
                       B6(:) = (B_face_zp(:,ii,jj,kk,lb)+B0_face_zp(:,ii,jj,kk,lb))*B_sc
                       
                       ele_nd =(G0_ion_prim(1,ii,jj,kk,lb)*rho_sc)/(M_p*amu)
                       rot_B(:) =  (/B4(3)-B3(3)+B5(2)-B6(2), B1(3)-B2(3)+B6(1)-B5(1), B2(2)-B1(2)+B3(1)-B4(1)/) 
                       Jd(:) = Rot_B(:)/(dx*l_sc*mu_0)
                       G0_ele_prim(2:4,ii,jj,kk,lb) = G0_ion_prim(2:4,ii,jj,kk,lb)-(Jd(:)/(ele_nd*q*a_sc))
                       J_prim(:,ii,jj,kk,lb) = Jd(:)/J_sc

                    end if

! Case of a surface cell -------------------------------------------------------!
                    if (loc_cell_flags(ii,jj,kk).gt.0) then

                       B1(:) = (B_face_xm(:,ii,jj,kk,lb)+B0_face_xm(:,ii,jj,kk,lb))*B_sc
                       B2(:) = (B_face_xp(:,ii,jj,kk,lb)+B0_face_xp(:,ii,jj,kk,lb))*B_sc
                       B3(:) = (B_face_ym(:,ii,jj,kk,lb)+B0_face_ym(:,ii,jj,kk,lb))*B_sc
                       B4(:) = (B_face_yp(:,ii,jj,kk,lb)+B0_face_yp(:,ii,jj,kk,lb))*B_sc
                       B5(:) = (B_face_zm(:,ii,jj,kk,lb)+B0_face_zm(:,ii,jj,kk,lb))*B_sc
                       B6(:) = (B_face_zp(:,ii,jj,kk,lb)+B0_face_zp(:,ii,jj,kk,lb))*B_sc
                       B7(:) = (B_face_ei(:,ii,jj,kk,lb)+B0_face_ei(:,ii,jj,kk,lb))*B_sc

                       rot_B(:) = 0.
                       do face=1,7                                              ! loop over cell's faces
                          if (surf_face_data(3,face,ii,jj,kk,bflags(1,lb)).ne.0) then
                             select case(face)
                             case(1)                                            ! compute flux on X-
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                rot_B(2) = rot_B(2)+(B1(3)*surface)
                                rot_B(3) = rot_B(3)-(B1(2)*surface)
                             case(2)                                            ! compute flux on X+
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                rot_B(2) = rot_B(2)-(B2(3)*surface)
                                rot_B(3) = rot_B(3)+(B2(2)*surface)
                             case(3)                                            ! compute flux on Y-
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                rot_B(1) = rot_B(1)-(B3(3)*surface)
                                rot_B(3) = rot_B(3)+(B3(1)*surface)
                            case(4)                                            ! compute flux on Y+
                               surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                               rot_B(1) = rot_B(1)+(B4(3)*surface)
                               rot_B(3) = rot_B(3)-(B4(1)*surface)
                            case(5)                                            ! compute flux on Z-
                               surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                               rot_B(1) = rot_B(1)+(B5(2)*surface)
                               rot_B(2) = rot_B(2)-(B5(1)*surface)
                             case(6)                                            ! compute flux on Z+ 
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                rot_B(1) = rot_B(1)-(B6(2)*surface)
                                rot_B(2) = rot_B(2)+(B6(1)*surface)
                             case(7)                                            ! compute flux on Cut face
                                nx = surf_face_data(4,face,ii,jj,kk,bflags(1,lb))
                                ny = surf_face_data(5,face,ii,jj,kk,bflags(1,lb))
                                nz = surf_face_data(6,face,ii,jj,kk,bflags(1,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                rot_B(1) = rot_B(1)+((ny*B7(3)-nz*B7(2))*surface)
                                rot_B(2) = rot_B(2)+((nz*B7(1)-nx*B7(3))*surface)
                                rot_B(3) = rot_B(3)+((nx*B7(2)-ny*B7(1))*surface)
                             end select
                          end if
                       end do

                       volume = surf_cell_data(1,ii,jj,kk,bflags(1,lb))*(dx_scale**3.)
                       ele_nd =(G0_ion_prim(1,ii,jj,kk,lb)*rho_sc)/(M_p*amu)
                       Jd(:)= Rot_B(:)/(volume*l_sc*mu_0)
                       J_face_xp(:,ii,jj,kk,lb) = Jd(:)
                       J_face_xm(:,ii,jj,kk,lb) = Jd(:)
                       J_face_yp(:,ii,jj,kk,lb) = Jd(:)
                       J_face_ym(:,ii,jj,kk,lb) = Jd(:)
                       J_face_zp(:,ii,jj,kk,lb) = Jd(:)
                       J_face_zm(:,ii,jj,kk,lb) = Jd(:)
                       J_face_ei(:,ii,jj,kk,lb) = Jd(:)
                       J_face_eo(:,ii,jj,kk,lb) = Jd(:)
                       G0_ele_prim(2:4,ii,jj,kk,lb) = G0_ion_prim(2:4,ii,jj,kk,lb)-(Jd(:)/(ele_nd*q*a_sc))
                       J_prim(:,ii,jj,kk,lb) = Jd(:)/J_sc

                    end if
                
                 end do
              end do
           end do

        else

           J_face_xp(:,:,:,:,lb) = 0.
           J_face_xm(:,:,:,:,lb) = 0.
           J_face_yp(:,:,:,:,lb) = 0.
           J_face_ym(:,:,:,:,lb) = 0.
           J_face_zp(:,:,:,:,lb) = 0.
           J_face_zm(:,:,:,:,lb) = 0.
           J_face_ei(:,:,:,:,lb) = 0.
           J_face_eo(:,:,:,:,lb) = 0.
           G0_ele_prim(2:4,:,:,:,lb) = G0_ion_prim(2:4,:,:,:,lb)
           J_prim(:,:,:,:,lb) = 0.

        end if

        G0_ele_prim(1,:,:,:,lb) = (G0_ion_prim(1,:,:,:,lb)*rho_sc)/(M_p*amu*nd_sc)
        G0_ele_cons(1,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)
        G0_ele_cons(2,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)*G0_ele_prim(2,:,:,:,lb)
        G0_ele_cons(3,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)*G0_ele_prim(3,:,:,:,lb)
        G0_ele_cons(4,:,:,:,lb) = G0_ele_prim(1,:,:,:,lb)*G0_ele_prim(4,:,:,:,lb)
        G0_ele_cons(5,:,:,:,lb) = G0_ele_prim(5,:,:,:,lb)
        J_save(1:3,1:nxb,1:nyb,1:nzb,lb) = J_prim(1:3,nguard+1:nxb+nguard,nguard+1:nyb+nguard,nguard+1:nzb+nguard,lb)

        return

        end subroutine current_face_surface

!-------------------------------------------------------------------------------!   
! Subroutine field_eignevalues                                                  !
!-------------------------------------------------------------------------------! 

        subroutine field_eigenvalues(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep

        call field_eigenvalues_x(lb)
        call field_eigenvalues_y(lb) 
        call field_eigenvalues_z(lb) 
        call field_eigenvalues_e(lb) 

        end subroutine field_eigenvalues

!-------------------------------------------------------------------------------!   
! Subroutine field_residual_normal                                              !
!-------------------------------------------------------------------------------! 

        subroutine field_residual_normal(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k, ii, jj, kk
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xm, flx_ym, flx_zm
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xp, flx_yp, flx_zp 
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_tx, flx_ty, flx_tz
        real, dimension(1:3,1:nxb,1:nyb,1:nzb) :: flx1, flx2, flx3, flx4, flx5, flx6
        real :: dx, lambda_x, lambda_y, lambda_z,dt_cell, dt_v
        real :: ion_vx, ion_vy, ion_vz, Bx1, Bx2, By3, By4, Bz5, Bz6
        real :: grad_B
        real, dimension(3) :: vec_B, div_term

! Compute cell size -------- ---------------------------------------------------!
        dx = (bsize(1,lb)*length_unit/(l_sc*nxb))

! Compute face fluxes ----------------------------------------------------------!
        call field_flux_function(lb,B_face_xm(:,:,:,:,lb), ele_face_xm(:,:,:,:,lb), ion_face_xm(:,:,:,:,lb), B0_face_xm(:,:,:,:,lb), flx_xm, flx_ty, flx_tz)
        call field_flux_function(lb,B_face_xp(:,:,:,:,lb), ele_face_xp(:,:,:,:,lb), ion_face_xp(:,:,:,:,lb), B0_face_xp(:,:,:,:,lb), flx_xp, flx_ty, flx_tz)
        call field_flux_function(lb,B_face_ym(:,:,:,:,lb), ele_face_ym(:,:,:,:,lb), ion_face_ym(:,:,:,:,lb), B0_face_ym(:,:,:,:,lb), flx_tx, flx_ym, flx_tz)
        call field_flux_function(lb,B_face_yp(:,:,:,:,lb), ele_face_yp(:,:,:,:,lb), ion_face_yp(:,:,:,:,lb), B0_face_yp(:,:,:,:,lb), flx_tx, flx_yp, flx_tz)
        call field_flux_function(lb,B_face_zm(:,:,:,:,lb), ele_face_zm(:,:,:,:,lb), ion_face_zm(:,:,:,:,lb), B0_face_zm(:,:,:,:,lb), flx_tx, flx_ty, flx_zm)
        call field_flux_function(lb,B_face_zp(:,:,:,:,lb), ele_face_zp(:,:,:,:,lb), ion_face_zp(:,:,:,:,lb), B0_face_zp(:,:,:,:,lb), flx_tx, flx_ty, flx_zp)

! Compute the total flux -------------------------------------------------------!
        do k=1,nzb                                                              ! loop over block's cells
           kk = k+nguard
           do j=1,nyb
              jj = j+nguard
              do i=1,nxb
                 ii= i+nguard

                 lambda_x = max(B_lambda_xp_av(ii,jj,kk,lb),B_lambda_xm_av(ii,jj,kk,lb))
                 lambda_y = max(B_lambda_yp_av(ii,jj,kk,lb),B_lambda_ym_av(ii,jj,kk,lb))
                 lambda_z = max(B_lambda_zp_av(ii,jj,kk,lb),B_lambda_zm_av(ii,jj,kk,lb))

                 flx1(:,i,j,k)=0.5*(flx_xp(:,ii-1,jj,kk)+flx_xm(:,ii,jj,kk))- & ! compute flux on X-
                               B_art*lambda_x*(B_face_xm(:,ii,jj,kk,lb)-B_face_xp(:,ii-1,jj,kk,lb))
                 flx2(:,i,j,k)=0.5*(flx_xp(:,ii,jj,kk)+flx_xm(:,ii+1,jj,kk))- & ! compute flux on X+
                               B_art*lambda_x*(B_face_xm(:,ii+1,jj,kk,lb)-B_face_xp(:,ii,jj,kk,lb))
                 flx3(:,i,j,k)=0.5*(flx_yp(:,ii,jj-1,kk)+flx_ym(:,ii,jj,kk))- & ! compute flux on Y-
                               B_art*lambda_y*(B_face_ym(:,ii,jj,kk,lb)-B_face_yp(:,ii,jj-1,kk,lb))
                 flx4(:,i,j,k)=0.5*(flx_yp(:,ii,jj,kk)+flx_ym(:,ii,jj+1,kk))- & ! compute flux on Y+
                               B_art*lambda_y*(B_face_ym(:,ii,jj+1,kk,lb)-B_face_yp(:,ii,jj,kk,lb))
                 flx5(:,i,j,k)=0.5*(flx_zp(:,ii,jj,kk-1)+flx_zm(:,ii,jj,kk))- & ! compute flux on Z-
                               B_art*lambda_z*(B_face_zm(:,ii,jj,kk,lb)-B_face_zp(:,ii,jj,kk-1,lb))
                 flx6(:,i,j,k)=0.5*(flx_zp(:,ii,jj,kk)+flx_zm(:,ii,jj,kk+1))- & ! compute flux on Z+
                               B_art*lambda_z*(B_face_zm(:,ii,jj,kk+1,lb)-B_face_zp(:,ii,jj,kk,lb))
                                                                           
! Load cell's values -----------------------------------------------------------!
                 ion_vx = G0_ion_prim(2,ii,jj,kk,lb)
                 ion_vy = G0_ion_prim(3,ii,jj,kk,lb)
                 ion_vz = G0_ion_prim(4,ii,jj,kk,lb)

! Compute the magnetic and velocity field gradient -----------------------------!
                 Bx1 = B_face_xm_av(1,ii,jj,kk,lb)
                 Bx2 = B_face_xp_av(1,ii,jj,kk,lb)
                 By3 = B_face_ym_av(2,ii,jj,kk,lb)
                 By4 = B_face_yp_av(2,ii,jj,kk,lb)
                 Bz5 = B_face_zm_av(3,ii,jj,kk,lb)
                 Bz6 = B_face_zp_av(3,ii,jj,kk,lb)
                 grad_B = (Bx2-Bx1+By4-By3+Bz6-Bz5)/dx
                 vec_B(:)=(/ion_vx, ion_vy, ion_vz/)
                 div_term(:) = vec_B(:)*grad_B
                 
! Compute ion residual ---------------------------------------------------------!
                 B_R(:,i,j,k,lb) = (flx6(:,i,j,k)-flx5(:,i,j,k)+flx4(:,i,j,k)-flx3(:,i,j,k)+flx2(:,i,j,k)-flx1(:,i,j,k))/dx + div_term(:)

                 if (init_sstep) then
                    
                    dt_v = CFL/((2.*(lambda_x+lambda_y+lambda_z)/dx)+(1./cell_time_stif(3,i,j,k,lb))+1.E-50)
                    dt_cell = min(dt_v,maxtimestep)
                    cell_time(3,i,j,k,lb) = dt_cell                             ! store cell's time step
                    local_time(3,lb) = min(dt_cell, local_time(3,lb))           ! compute block's field time step 
                    local_time(4,lb) = min(dt_cell, local_time(4,lb))           ! compute block's absolute time step 

                 end if

              end do
           end do
        end do

        return

        end subroutine field_residual_normal

!-------------------------------------------------------------------------------!   
! Subroutine field_adv_normal                                                   !
!-------------------------------------------------------------------------------! 

        subroutine field_adv_normal(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k, ii, jj, kk
        real :: dt
        real, dimension(1:nxb,1:nyb,1:nzb) :: res

        res(:,:,:) = 0.

! Advance the solution ---------------------------------------------------------!
        do k = nguard+1,nzb+nguard                                              ! loop over the block's cells
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
                                                                                ! compute next solution
                 G1_B_prim(:,i,j,k,lb)=G0_B_prim(:,i,j,k,lb)-(dt*B_R(:,ii,jj,kk,lb)) 
                 
                 res(ii,jj,kk)=sqrt(B_R(1,ii,jj,kk,lb)**2.+B_R(2,ii,jj,kk,lb)**2.+B_R(3,ii,jj,kk,lb)**2.)
      
              end do
           end do
        end do
                       
! Compute block's local residual -----------------------------------------------!
        local_residu(3,lb) = maxval(res)
     
        return

        end subroutine field_adv_normal

!-------------------------------------------------------------------------------!   
! Subroutine field_residual_surface                                             !
!-------------------------------------------------------------------------------! 

        subroutine field_residual_surface(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep

        integer :: i, j, k, ii, jj, kk, face
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xm, flx_ym, flx_zm
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_xp, flx_yp, flx_zp 
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_tx, flx_ty, flx_tz
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_eo_xm, flx_eo_ym, flx_eo_zm
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: flx_ei_xp, flx_ei_yp, flx_ei_zp
        real, dimension(1:3,1:nxb,1:nyb,1:nzb) :: flx1, flx2, flx3, flx4, flx5, flx6, flx7, flxa, flxb, flxc, flx, flx_tot
        real :: charac, lambda_x_e, lambda_y_e, lambda_z_e, lambda_x, lambda_y, lambda_z
        real :: ion_vx, ion_vy, ion_vz, Bx1, Bx2, By3, By4, Bz5, Bz6, surface, nx, ny, nz, volume
        real :: Vx1, Vx2, Vy3, Vy4, Vz5, Vz6, grad_B
        real, dimension(3) :: vec_B, div_term
        real :: dt_cell, dt_v, dx_scale, dx
 
! Compute cell size -------- ---------------------------------------------------!
        dx = (bsize(1,lb)*length_unit/(l_sc*nxb))
        dx_scale = length_unit/l_sc

! Compute face fluxes ----------------------------------------------------------!
        call field_flux_function(lb,B_face_xm(:,:,:,:,lb), ele_face_xm(:,:,:,:,lb), ion_face_xm(:,:,:,:,lb), B0_face_xm(:,:,:,:,lb), flx_xm, flx_ty, flx_tz)
        call field_flux_function(lb,B_face_xp(:,:,:,:,lb), ele_face_xp(:,:,:,:,lb), ion_face_xp(:,:,:,:,lb), B0_face_xp(:,:,:,:,lb), flx_xp, flx_ty, flx_tz)
        call field_flux_function(lb,B_face_ym(:,:,:,:,lb), ele_face_ym(:,:,:,:,lb), ion_face_ym(:,:,:,:,lb), B0_face_ym(:,:,:,:,lb), flx_tx, flx_ym, flx_tz)
        call field_flux_function(lb,B_face_yp(:,:,:,:,lb), ele_face_yp(:,:,:,:,lb), ion_face_yp(:,:,:,:,lb), B0_face_yp(:,:,:,:,lb), flx_tx, flx_yp, flx_tz)
        call field_flux_function(lb,B_face_zm(:,:,:,:,lb), ele_face_zm(:,:,:,:,lb), ion_face_zm(:,:,:,:,lb), B0_face_zm(:,:,:,:,lb), flx_tx, flx_ty, flx_zm)
        call field_flux_function(lb,B_face_zp(:,:,:,:,lb), ele_face_zp(:,:,:,:,lb), ion_face_zp(:,:,:,:,lb), B0_face_zp(:,:,:,:,lb), flx_tx, flx_ty, flx_zp)
        call field_flux_function(lb,B_face_eo(:,:,:,:,lb), ele_face_eo(:,:,:,:,lb), ion_face_eo(:,:,:,:,lb), B0_face_eo(:,:,:,:,lb), flx_eo_xm, flx_eo_ym, flx_eo_zm)
        call field_flux_function(lb,B_face_ei(:,:,:,:,lb), ele_face_ei(:,:,:,:,lb), ion_face_ei(:,:,:,:,lb), B0_face_ei(:,:,:,:,lb), flx_ei_xp, flx_ei_yp, flx_ei_zp)

! Compute the flux functions ---------------------------------------------------!
        do k=1,nzb                                                              ! loop over block's cells
           kk = k+nguard
           do j=1,nyb
              jj = j+nguard
              do i=1,nxb
                 ii= i+nguard

! Case of a non surface cell ---------------------------------------------------!
                 if (loc_cell_flags(ii,jj,kk).eq.0) then

                    lambda_x = max(B_lambda_xp_av(ii,jj,kk,lb),B_lambda_xm_av(ii,jj,kk,lb))
                    lambda_y = max(B_lambda_yp_av(ii,jj,kk,lb),B_lambda_ym_av(ii,jj,kk,lb))
                    lambda_z = max(B_lambda_zp_av(ii,jj,kk,lb),B_lambda_zm_av(ii,jj,kk,lb))

                    flx1(:,i,j,k)=0.5*(flx_xp(:,ii-1,jj,kk)+flx_xm(:,ii,jj,kk))-& ! compute flux on X-
                         B_art*lambda_x*(B_face_xm(:,ii,jj,kk,lb)-B_face_xp(:,ii-1,jj,kk,lb))
                    flx2(:,i,j,k)=0.5*(flx_xp(:,ii,jj,kk)+flx_xm(:,ii+1,jj,kk))-& ! compute flux on X+
                         B_art*lambda_x*(B_face_xm(:,ii+1,jj,kk,lb)-B_face_xp(:,ii,jj,kk,lb))
                    flx3(:,i,j,k)=0.5*(flx_yp(:,ii,jj-1,kk)+flx_ym(:,ii,jj,kk))-& ! compute flux on Y-
                         B_art*lambda_y*(B_face_ym(:,ii,jj,kk,lb)-B_face_yp(:,ii,jj-1,kk,lb))
                    flx4(:,i,j,k)=0.5*(flx_yp(:,ii,jj,kk)+flx_ym(:,ii,jj+1,kk))-& ! compute flux on Y+
                         B_art*lambda_y*(B_face_ym(:,ii,jj+1,kk,lb)-B_face_yp(:,ii,jj,kk,lb))
                    flx5(:,i,j,k)=0.5*(flx_zp(:,ii,jj,kk-1)+flx_zm(:,ii,jj,kk))-& ! compute flux on Z-
                         B_art*lambda_z*(B_face_zm(:,ii,jj,kk,lb)-B_face_zp(:,ii,jj,kk-1,lb))
                    flx6(:,i,j,k)=0.5*(flx_zp(:,ii,jj,kk)+flx_zm(:,ii,jj,kk+1))-& ! compute flux on Z+
                         B_art*lambda_z*(B_face_zm(:,ii,jj,kk+1,lb)-B_face_zp(:,ii,jj,kk,lb)) 
              
! Load cell's values -----------------------------------------------------------!                              
                    ion_vx = G0_ion_prim(2,ii,jj,kk,lb)
                    ion_vy = G0_ion_prim(3,ii,jj,kk,lb)
                    ion_vz = G0_ion_prim(4,ii,jj,kk,lb)

! Compute the magnetic and velocity field gradient -----------------------------!
                    Bx1 = B_face_xm_av(1,ii,jj,kk,lb)
                    Bx2 = B_face_xp_av(1,ii,jj,kk,lb)
                    By3 = B_face_ym_av(2,ii,jj,kk,lb)
                    By4 = B_face_yp_av(2,ii,jj,kk,lb)
                    Bz5 = B_face_zm_av(3,ii,jj,kk,lb)
                    Bz6 = B_face_zp_av(3,ii,jj,kk,lb)
                    grad_B = (Bx2-Bx1+By4-By3+Bz6-Bz5)/dx
                    vec_B(:)=(/ion_vx, ion_vy, ion_vz/)
                    div_term(:) = vec_B(:)*grad_B

! Compute ion residual ---------------------------------------------------------!
                    B_R(:,i,j,k,lb) = (flx6(:,i,j,k)-flx5(:,i,j,k)+flx4(:,i,j,k)-flx3(:,i,j,k)+flx2(:,i,j,k)-flx1(:,i,j,k))/dx + div_term(:)

                    if (init_sstep) then

                       dt_v = CFL/((2.*(lambda_x+lambda_y+lambda_z)/dx)+(1./cell_time_stif(3,i,j,k,lb))+1.E-50)
                       dt_cell = min(dt_v,maxtimestep)
                       cell_time(3,i,j,k,lb) = dt_cell                          ! store cell's time step
                       local_time(3,lb) = min(dt_cell, local_time(3,lb))        ! compute block's field time step 
                       local_time(4,lb) = min(dt_cell, local_time(4,lb))        ! compute block's absolute time step 

                    end if

                 end if

! Case of a surface cell -------------------------------------------------------!
                 if (loc_cell_flags(ii,jj,kk).gt.0) then
                                                                                ! compute average eignvalues 
                    lambda_x = 0.
                    lambda_y = 0.
                    lambda_z = 0.
                    if (surf_face_data(3,1,ii,jj,kk,bflags(1,lb)).ne.0) lambda_x = max(lambda_x,B_lambda_xm_av(ii,jj,kk,lb))
                    if (surf_face_data(3,2,ii,jj,kk,bflags(1,lb)).ne.0) lambda_x = max(lambda_x,B_lambda_xp_av(ii,jj,kk,lb))
                    if (surf_face_data(3,3,ii,jj,kk,bflags(1,lb)).ne.0) lambda_y = max(lambda_y,B_lambda_ym_av(ii,jj,kk,lb))
                    if (surf_face_data(3,4,ii,jj,kk,bflags(1,lb)).ne.0) lambda_y = max(lambda_y,B_lambda_yp_av(ii,jj,kk,lb))
                    if (surf_face_data(3,5,ii,jj,kk,bflags(1,lb)).ne.0) lambda_z = max(lambda_z,B_lambda_zm_av(ii,jj,kk,lb))
                    if (surf_face_data(3,6,ii,jj,kk,bflags(1,lb)).ne.0) lambda_z = max(lambda_z,B_lambda_zp_av(ii,jj,kk,lb))
                    if (surf_face_data(3,7,ii,jj,kk,bflags(1,lb)).ne.0) then
                       lambda_x = max(lambda_x,B_lambda_ex_av(ii,jj,kk,lb))
                       lambda_y = max(lambda_y,B_lambda_ey_av(ii,jj,kk,lb))
                       lambda_z = max(lambda_z,B_lambda_ez_av(ii,jj,kk,lb))
                    end if

                    charac = 0.
                    flx_tot(:,:,:,:) = 0.
                    grad_B = 0.

                    do face=1,7                                                 ! loop over cell's faces
                       if (surf_face_data(3,face,ii,jj,kk,bflags(1,lb)).ne.0) then
                          select case(face)
                             case(1)                                            ! compute flux on X-
                                flx1(:,i,j,k)=0.5*(flx_xp(:,ii-1,jj,kk)+flx_xm(:,ii,jj,kk))- &
                                     B_art*lambda_x*(B_face_xm(:,ii,jj,kk,lb)-B_face_xp(:,ii-1,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)-(flx1(:,i,j,k)*surface)
                                charac = charac+lambda_x*surface
                                grad_B =grad_B-(B_face_xm_av(1,ii,jj,kk,lb)*surface)
                             case(2)                                            ! compute flux on X+
                                flx2(:,i,j,k)=0.5*(flx_xp(:,ii,jj,kk)+flx_xm(:,ii+1,jj,kk))- &
                                     B_art*lambda_x*(B_face_xm(:,ii+1,jj,kk,lb)-B_face_xp(:,ii,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)+(flx2(:,i,j,k)*surface)
                                charac = charac+lambda_x*surface
                                grad_B =grad_B+(B_face_xp_av(1,ii,jj,kk,lb)*surface)
                             case(3)                                            ! compute flux on Y-
                                flx3(:,i,j,k)=0.5*(flx_yp(:,ii,jj-1,kk)+flx_ym(:,ii,jj,kk))- &
                                     B_art*lambda_y*(B_face_ym(:,ii,jj,kk,lb)-B_face_yp(:,ii,jj-1,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)-(flx3(:,i,j,k)*surface)
                                charac = charac+lambda_y*surface
                                grad_B =grad_B-(B_face_ym_av(2,ii,jj,kk,lb)*surface)
                             case(4)                                            ! compute flux on Y+
                                flx4(:,i,j,k)=0.5*(flx_yp(:,ii,jj,kk)+flx_ym(:,ii,jj+1,kk))- &
                                     B_art*lambda_y*(B_face_ym(:,ii,jj+1,kk,lb)-B_face_yp(:,ii,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)+(flx4(:,i,j,k)*surface)
                                charac = charac+lambda_y*surface
                                grad_B =grad_B+(B_face_yp_av(2,ii,jj,kk,lb)*surface)
                             case(5)                                            ! compute flux on Z-
                                flx5(:,i,j,k)=0.5*(flx_zp(:,ii,jj,kk-1)+flx_zm(:,ii,jj,kk))- &
                                     B_art*lambda_z*(B_face_zm(:,ii,jj,kk,lb)-B_face_zp(:,ii,jj,kk-1,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)-(flx5(:,i,j,k)*surface)
                                charac = charac+lambda_z*surface
                                grad_B =grad_B-(B_face_zm_av(3,ii,jj,kk,lb)*surface)
                             case(6)                                            ! compute flux on Z+ 
                                flx6(:,i,j,k)=0.5*(flx_zp(:,ii,jj,kk)+flx_zm(:,ii,jj,kk+1))- &
                                     B_art*lambda_z*(B_face_zm(:,ii,jj,kk+1,lb)-B_face_zp(:,ii,jj,kk,lb))
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)+(flx6(:,i,j,k)*surface)
                                charac = charac+lambda_z*surface
                                grad_B =grad_B+(B_face_zp_av(3,ii,jj,kk,lb)*surface)
                             case(7)                                            ! compute flux on Cut face
                                nx = surf_face_data(4,face,ii,jj,kk,bflags(1,lb))
                                ny = surf_face_data(5,face,ii,jj,kk,bflags(1,lb))
                                nz = surf_face_data(6,face,ii,jj,kk,bflags(1,lb))
                                flxa(:,i,j,k)=0.5*(flx_ei_xp(:,ii,jj,kk)+flx_eo_xm(:,ii,jj,kk))- &
                                     sign(1.,nx)*B_art*lambda_x*(B_face_eo(:,ii,jj,kk,lb)-B_face_ei(:,ii,jj,kk,lb))
                                flxb(:,i,j,k)=0.5*(flx_ei_yp(:,ii,jj,kk)+flx_eo_ym(:,ii,jj,kk))- &
                                     sign(1.,ny)*B_art*lambda_y*(B_face_eo(:,ii,jj,kk,lb)-B_face_ei(:,ii,jj,kk,lb))
                                flxc(:,i,j,k)=0.5*(flx_ei_zp(:,ii,jj,kk)+flx_eo_zm(:,ii,jj,kk))- &
                                     sign(1.,nz)*B_art*lambda_z*(B_face_eo(:,ii,jj,kk,lb)-B_face_ei(:,ii,jj,kk,lb))
                                flx7(:,i,j,k) = (flxa(:,i,j,k)*nx+flxb(:,i,j,k)*ny+flxc(:,i,j,k)*nz)
                                surface = surf_face_data(7,face,ii,jj,kk,bflags(1,lb))*dx_scale*dx_scale
                                flx_tot(:,i,j,k)=flx_tot(:,i,j,k)+(flx7(:,i,j,k)*surface)
                                charac = charac+(lambda_x+lambda_y+lambda_z)*surface
                                grad_B =grad_B+(B_face_ex_av(1,ii,jj,kk,lb)*nx+B_face_ex_av(2,ii,jj,kk,lb)*ny &
                                               +B_face_ex_av(3,ii,jj,kk,lb)*nz)*surface
                          end select
                       end if
                    end do
                                                                                ! compute cell's volume
                    volume = surf_cell_data(1,ii,jj,kk,bflags(1,lb))*(dx_scale**3.)

! Load cell's values -----------------------------------------------------------!
                    ion_vx = G0_ion_prim(2,ii,jj,kk,lb)
                    ion_vy = G0_ion_prim(3,ii,jj,kk,lb)
                    ion_vz = G0_ion_prim(4,ii,jj,kk,lb)
                    
! Compute the magnetic and velocity field --------------------------------------!
                    grad_B=grad_B/volume
                    vec_B(:)=(/ion_vx, ion_vy, ion_vz/)
                    div_term(:) = vec_B(:)*grad_B

! Compute ion residual ---------------------------------------------------------!
                    B_R(:,i,j,k,lb)=(flx_tot(:,i,j,k)/volume)+div_term(:)

                    if (init_sstep) then

                       dt_v = (CFL*volume)/(charac+(volume/cell_time_stif(3,i,j,k,lb))+1.e-50)
                       dt_cell = min(dt_v,maxtimestep)
                       cell_time(3,i,j,k,lb) = dt_cell                          ! store cell's time step
                       local_time(3,lb) = min(dt_cell, local_time(3,lb))        ! compute block's field time step 

                    end if
   
                 end if
                 
              end do
           end do
        end do

        return
        
        end subroutine field_residual_surface 

!-------------------------------------------------------------------------------!   
! Subroutine field_adv_surface                                                  !
!-------------------------------------------------------------------------------! 

        subroutine field_adv_surface(mype, lb, istep)

        implicit none 

        integer, intent(in) ::lb, mype, istep
        integer :: i, j, k, ii, jj, kk
        real :: dt
        real, dimension(1:nxb,1:nyb,1:nzb) :: res

        res(:,:,:) = 0.
        G1_B_prim(:,:,:,:,lb) = 2.
! Advance the solution ---------------------------------------------------------!
        do k = nguard+1,nzb+nguard                                              ! loop over the block's cells
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

                    G1_B_prim(:,i,j,k,lb)=G0_B_prim(:,i,j,k,lb)-(dt*B_R(:,ii,jj,kk,lb))
                    
                    res(ii,jj,kk)=sqrt(B_R(1,ii,jj,kk,lb)**2.+B_R(2,ii,jj,kk,lb)**2.+B_R(3,ii,jj,kk,lb)**2.)
 
                 end if

              end do
           end do
        end do

! Compute block's local residual -----------------------------------------------!
        local_residu(3,lb) = maxval(res)

        return

        end subroutine field_adv_surface
          
!-------------------------------------------------------------------------------!   
! Subroutine field_interpol_to_faces_normal                                     !
!-------------------------------------------------------------------------------! 

        subroutine field_interpol_to_faces_normal(lb, istep, vec, fxm, fxp, fym, fyp, fzm, fzp, ord)

        implicit none 
        
        integer, intent(in) :: lb, istep, ord
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(inout) :: vec
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: fxm, fxp, fym, fyp, fzm, fzp
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: fei
        integer :: i, j, k

! Reconstruct cell edges -------------------------------------------------------!
        do k=2,nzb+2*nguard-1
           do j=2,nyb+2*nguard-1
              do i=2,nxb+2*nguard-1
                 call field_reconstruct_edge(lb, vec(:,i,j,k), vec(:,i-1,j,k), vec(:,i+1,j,k), vec(:,i,j-1,k), vec(:,i,j+1,k), &
                 vec(:,i,j,k-1), vec(:,i,j,k+1), vec(:,i-1,j-1,k), vec(:,i+1,j-1,k), vec(:,i+1,j+1,k), vec(:,i-1,j+1,k), &
                 vec(:,i-1,j-1,k-1), vec(:,i+1,j-1,k-1), vec(:,i+1,j+1,k-1), vec(:,i-1,j+1,k-1), vec(:,i-1,j-1,k+1), vec(:,i+1,j-1,k+1), &
                 vec(:,i+1,j+1,k+1), vec(:,i-1,j+1,k+1), vec(:,i-1,j,k-1), vec(:,i+1,j,k-1), vec(:,i,j-1,k-1), vec(:,i,j+1,k-1), &
                 vec(:,i-1,j,k+1), vec(:,i+1,j,k+1), vec(:,i,j-1,k+1), vec(:,i,j+1,k+1), fxm(:,i,j,k), fxp(:,i,j,k), fym(:,i,j,k), &
                 fyp(:,i,j,k), fzm(:,i,j,k), fzp(:,i,j,k), fei(:,i,j,k), ord, dxcell(i,j,k,lb))

              end do
           end do
        end do

        return
        end subroutine field_interpol_to_faces_normal

!-------------------------------------------------------------------------------!   
! Subroutine field_interpol_to_faces_surface                                    !
!-------------------------------------------------------------------------------! 

        subroutine field_interpol_to_faces_surface(lb, istep, vec, fxm, fxp, fym, fyp, fzm, fzp, fei, feo, ord)

        implicit none 
        
        integer, intent(in) :: lb, istep, ord
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(inout) :: vec
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: fxm, fxp, fym, fyp, fzm, fzp, fei, feo
        integer :: i, j, k, ii, jj, kk, face, nuorder
        real :: xc, yc, zc
 
! Reconstruct cell edges -------------------------------------------------------!
        do k=2,nzb+2*nguard-1
           do j=2,nyb+2*nguard-1
              do i=2,nxb+2*nguard-1

                 if (loc_cell_flags(i,j,k).ne.-1) then                          ! check that it is not a nucleus cell

                    if(loc_cell_flags(i,j,k).eq.0) then                         ! case of non surface cell
                       nuorder = ord
                       call field_reconstruct_edge(lb, vec(:,i,j,k), vec(:,i-1,j,k), vec(:,i+1,j,k), vec(:,i,j-1,k), vec(:,i,j+1,k), &
                       vec(:,i,j,k-1), vec(:,i,j,k+1), vec(:,i-1,j-1,k), vec(:,i+1,j-1,k), vec(:,i+1,j+1,k), vec(:,i-1,j+1,k), &
                       vec(:,i-1,j-1,k-1), vec(:,i+1,j-1,k-1), vec(:,i+1,j+1,k-1), vec(:,i-1,j+1,k-1), vec(:,i-1,j-1,k+1), &
                       vec(:,i+1,j-1,k+1), vec(:,i+1,j+1,k+1), vec(:,i-1,j+1,k+1), vec(:,i-1,j,k-1), vec(:,i+1,j,k-1), vec(:,i,j-1,k-1), &
                       vec(:,i,j+1,k-1),  vec(:,i-1,j,k+1), vec(:,i+1,j,k+1), vec(:,i,j-1,k+1), vec(:,i,j+1,k+1), fxm(:,i,j,k), &
                       fxp(:,i,j,k), fym(:,i,j,k), fyp(:,i,j,k), fzm(:,i,j,k),fzp(:,i,j,k), fei(:,i,j,k), nuorder, dxcell(i,j,k,lb))
                    end if

                    if (loc_cell_flags(i,j,k).gt.0) then                        ! case of a surface cell
                       nuorder = 2   
                       call field_reconstruct_edge2(lb, i, j, k,  vec(:,i,j,k), vec(:,i-1,j,k), vec(:,i+1,j,k), vec(:,i,j-1,k), &
                       vec(:,i,j+1,k), vec(:,i,j,k-1), vec(:,i,j,k+1), vec(:,i-1,j-1,k), vec(:,i+1,j-1,k), vec(:,i+1,j+1,k), &
                       vec(:,i-1,j+1,k), vec(:,i-1,j-1,k-1), vec(:,i+1,j-1,k-1), vec(:,i+1,j+1,k-1), vec(:,i-1,j+1,k-1), &
                       vec(:,i-1,j-1,k+1), vec(:,i+1,j-1,k+1), vec(:,i+1,j+1,k+1), vec(:,i-1,j+1,k+1), vec(:,i-1,j,k-1), &
                       vec(:,i+1,j,k-1), vec(:,i,j-1,k-1), vec(:,i,j+1,k-1), vec(:,i-1,j,k+1), vec(:,i+1,j,k+1), vec(:,i,j-1,k+1), &
                       vec(:,i,j+1,k+1), fxm(:,i,j,k), fxp(:,i,j,k), fym(:,i,j,k), fyp(:,i,j,k), fzm(:,i,j,k), fzp(:,i,j,k), &
                       fei(:,i,j,k), nuorder)

                       if (surf_face_data(3,7,i,j,k,bflags(1,lb)).ne.0) then    ! apply boundary condition
                          xc = surf_face_data(8,7,i,j,k,bflags(1,lb))*length_unit
                          yc = surf_face_data(9,7,i,j,k,bflags(1,lb))*length_unit
                          zc = surf_face_data(10,7,i,j,k,bflags(1,lb))*length_unit
                          call field_set_surf_bc(lb, xc, yc, zc, fei(:,i,j,k), feo(:,i,j,k))
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
        end subroutine  field_interpol_to_faces_surface

!-------------------------------------------------------------------------------!   
! Subroutine field_set_surf_bc                                                  !
!-------------------------------------------------------------------------------!  
        subroutine field_set_surf_bc(lb, xc, yc, zc, fi, fo)

        implicit none 
        integer, intent(in) :: lb
        real, intent(in) :: xc, yc, zc
        real, dimension(1:3), intent(inout) :: fi, fo

        fo(:) = 0.

        return
        end subroutine field_set_surf_bc

!-------------------------------------------------------------------------------!   
! Subroutine field_eigenvalues_x                                                !
!-------------------------------------------------------------------------------!

        subroutine field_eigenvalues_x(lb)

        implicit none 
        integer, intent(in) :: lb
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: Bx, By, Bz, Jx, Jy, Jz, cf, n, cw

        Bx(:,:,:) = B_face_xm_av(1,:,:,:,lb)+B0_face_xm(1,:,:,:,lb)
        By(:,:,:) = B_face_xm_av(2,:,:,:,lb)+B0_face_xm(2,:,:,:,lb)
        Bz(:,:,:) = B_face_xm_av(3,:,:,:,lb)+B0_face_xm(3,:,:,:,lb)
        n(:,:,:)= ele_face_xm_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if
                    
        B_lambda_xm_av(:,:,:,lb) = ion_lambda_xm_av(:,:,:,lb)+cw(:,:,:)

        Bx(:,:,:) = B_face_xp_av(1,:,:,:,lb)+B0_face_xp(1,:,:,:,lb)
        By(:,:,:) = B_face_xp_av(2,:,:,:,lb)+B0_face_xp(2,:,:,:,lb)
        Bz(:,:,:) = B_face_xp_av(3,:,:,:,lb)+B0_face_xp(3,:,:,:,lb)
        n(:,:,:)= ele_face_xp_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if

        B_lambda_xp_av(:,:,:,lb) = ion_lambda_xp_av(:,:,:,lb)+cw(:,:,:)

        return

        end subroutine field_eigenvalues_x

!-------------------------------------------------------------------------------!   
! Subroutine field_eigenvalues_y                                                !
!-------------------------------------------------------------------------------!

        subroutine field_eigenvalues_y(lb)

        implicit none
        integer, intent(in) :: lb
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: Bx, By, Bz, Jx, Jy, Jz,  n, cw
      
        Bx(:,:,:) = B_face_ym_av(1,:,:,:,lb)+B0_face_ym(1,:,:,:,lb)
        By(:,:,:) = B_face_ym_av(2,:,:,:,lb)+B0_face_ym(2,:,:,:,lb)
        Bz(:,:,:) = B_face_ym_av(3,:,:,:,lb)+B0_face_ym(3,:,:,:,lb)
        n(:,:,:)= ele_face_ym_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if

        B_lambda_ym_av(:,:,:,lb) = ion_lambda_ym_av(:,:,:,lb)+cw(:,:,:)
     
        Bx(:,:,:) = B_face_yp_av(1,:,:,:,lb)+B0_face_yp(1,:,:,:,lb)
        By(:,:,:) = B_face_yp_av(2,:,:,:,lb)+B0_face_yp(2,:,:,:,lb)
        Bz(:,:,:) = B_face_yp_av(3,:,:,:,lb)+B0_face_yp(3,:,:,:,lb)
        n(:,:,:)= ele_face_yp_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if

        B_lambda_yp_av(:,:,:,lb) = ion_lambda_yp_av(:,:,:,lb)+cw(:,:,:)

        return

        end subroutine field_eigenvalues_y

!-------------------------------------------------------------------------------!   
! Subroutine field_eigenvalues_z                                                !
!-------------------------------------------------------------------------------!

        subroutine field_eigenvalues_z(lb)

        implicit none
        integer, intent(in) :: lb
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: Bx, By, Bz, Jx, Jy, Jz, n, cw

        Bx(:,:,:) = B_face_zm_av(1,:,:,:,lb)+B0_face_zm(1,:,:,:,lb)
        By(:,:,:) = B_face_zm_av(2,:,:,:,lb)+B0_face_zm(2,:,:,:,lb)
        Bz(:,:,:) = B_face_zm_av(3,:,:,:,lb)+B0_face_zm(3,:,:,:,lb)
        n(:,:,:)= ele_face_zm_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if
        
        B_lambda_zm_av(:,:,:,lb) = ion_lambda_zm_av(:,:,:,lb)+cw(:,:,:)

        Bx(:,:,:) = B_face_zp_av(1,:,:,:,lb)+B0_face_zp(1,:,:,:,lb)
        By(:,:,:) = B_face_zp_av(2,:,:,:,lb)+B0_face_zp(2,:,:,:,lb)
        Bz(:,:,:) = B_face_zp_av(3,:,:,:,lb)+B0_face_zp(3,:,:,:,lb)
        n(:,:,:)= ele_face_zp_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if

        B_lambda_zp_av(:,:,:,lb) = ion_lambda_zp_av(:,:,:,lb)+cw(:,:,:)

        return

        end subroutine field_eigenvalues_z

!-------------------------------------------------------------------------------!   
! Subroutine field_eigenvalues_e                                                !
!-------------------------------------------------------------------------------!

        subroutine field_eigenvalues_e(lb)

        implicit none
        integer, intent(in) :: lb
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) ::Bx, By, Bz, Jx, Jy, Jz, n, cw
       
        Bx(:,:,:) = B_face_ex_av(1,:,:,:,lb)+B0_face_ei(1,:,:,:,lb)
        By(:,:,:) = B_face_ex_av(2,:,:,:,lb)+B0_face_ei(2,:,:,:,lb)
        Bz(:,:,:) = B_face_ex_av(3,:,:,:,lb)+B0_face_ei(3,:,:,:,lb)
        n(:,:,:)= ele_face_ex_av(1,:,:,:,lb)

        if (hall_mhd_flag==.true.) then
           cw(:,:,:) = 2.*(sqrt(Bx(:,:,:)**2.+By(:,:,:)**2.+Bz(:,:,:)**2.)*B_sc*pi)/(mu_0*q*n(:,:,:)*nd_sc*dxcell(:,:,:,lb)*l_sc*a_sc)
        else
           cw(:,:,:)=0.
        end if

        B_lambda_ex_av(:,:,:,lb) = ion_lambda_ex_av(:,:,:,lb)+cw(:,:,:)
        B_lambda_ey_av(:,:,:,lb) = ion_lambda_ey_av(:,:,:,lb)+cw(:,:,:)
        B_lambda_ez_av(:,:,:,lb) = ion_lambda_ez_av(:,:,:,lb)+cw(:,:,:)

        return

        end subroutine field_eigenvalues_e

!-------------------------------------------------------------------------------!   
! Subroutine field_adjust_storage                                               !
!-------------------------------------------------------------------------------!

        subroutine field_adjust_storage(lb, conserv)

        implicit none

        integer, intent(in) :: lb
        real, dimension(1:B_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: conserv
        integer :: i ,j, k

        do k=nguard+1,nzb+nguard
           do j=nguard+1,nyb+nguard
              do i=nguard+1,nxb+nguard 

                 if (ieee_is_nan(conserv(1,i,j,k)).or.ieee_is_nan(conserv(2,i,j,k)).or.ieee_is_nan(conserv(3,i,j,k))) then
                    exception_flag = 1
                    print*,'Exception detected Field_Module :', trans_mype, lb, i,j,k, conserv(:,i,j,k)
                 end if

              end do
           end do
        end do

        return
        end subroutine field_adjust_storage

!-------------------------------------------------------------------------------!   
! Subroutine field_flux_function                                                  !
!-------------------------------------------------------------------------------!

        subroutine field_flux_function(lb, B_prim, ele_prim, ion_prim, B0, flx_x, flx_y, flx_z)

        implicit none
        integer, intent(in) :: lb
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: B_prim
        real, dimension(1:ele_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: ele_prim
        real, dimension(1:ion_nvar,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: ion_prim
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(in) :: B0
        real, dimension(1:3,1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard), intent(out):: flx_x, flx_y, flx_z
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: Bx, By, Bz, B0x, B0y, B0z
        real, dimension(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard) :: vex, vey, vez

        Bx(:,:,:) = B_prim(1,:,:,:)
        By(:,:,:) = B_prim(2,:,:,:)
        Bz(:,:,:) = B_prim(3,:,:,:)

        B0x(:,:,:) = B0(1,:,:,:)
        B0y(:,:,:) = B0(2,:,:,:)
        B0z(:,:,:) = B0(3,:,:,:)

        vex(:,:,:) = ele_prim(2,:,:,:)
        vey(:,:,:) = ele_prim(3,:,:,:)
        vez(:,:,:) = ele_prim(4,:,:,:)

        flx_x(1,:,:,:) = 0.
        flx_y(1,:,:,:) = -(vex(:,:,:)*By(:,:,:))+(Bx(:,:,:)*vey(:,:,:))
        flx_z(1,:,:,:) = -(vex(:,:,:)*Bz(:,:,:))+(Bx(:,:,:)*vez(:,:,:))

        flx_x(2,:,:,:) = -(vey(:,:,:)*Bx(:,:,:))+(By(:,:,:)*vex(:,:,:))
        flx_y(2,:,:,:) = 0.
        flx_z(2,:,:,:) = -(vey(:,:,:)*Bz(:,:,:))+(By(:,:,:)*vez(:,:,:))

        flx_x(3,:,:,:) = -(vez(:,:,:)*Bx(:,:,:))+(Bz(:,:,:)*vex(:,:,:))
        flx_y(3,:,:,:) = -(vez(:,:,:)*By(:,:,:))+(Bz(:,:,:)*vey(:,:,:))
        flx_z(3,:,:,:) = 0.

! Compute the B0 contribution --------------------------------------------------!

        flx_y(1,:,:,:) = flx_y(1,:,:,:)-(vex(:,:,:)*B0y(:,:,:))+(B0x(:,:,:)*vey(:,:,:))
        flx_z(1,:,:,:) = flx_z(1,:,:,:)-(vex(:,:,:)*B0z(:,:,:))+(B0x(:,:,:)*vez(:,:,:))

        flx_x(2,:,:,:) = flx_x(2,:,:,:)-(vey(:,:,:)*B0x(:,:,:))+(B0y(:,:,:)*vex(:,:,:))
        flx_z(2,:,:,:) = flx_z(2,:,:,:)-(vey(:,:,:)*B0z(:,:,:))+(B0y(:,:,:)*vez(:,:,:))

        flx_x(3,:,:,:) = flx_x(3,:,:,:)-(vez(:,:,:)*B0x(:,:,:))+(B0z(:,:,:)*vex(:,:,:))
        flx_y(3,:,:,:) = flx_y(3,:,:,:)-(vez(:,:,:)*B0y(:,:,:))+(B0z(:,:,:)*vey(:,:,:))

        return
        end subroutine field_flux_function

!-------------------------------------------------------------------------------!   
! Subroutine field_reconstruct_edge                                             !
!-------------------------------------------------------------------------------!

        subroutine field_reconstruct_edge(lb, W_0, W_1, W_2, W_3, W_4, W_5, W_6, &
        C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12, E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, &
        W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7, ord, dx)

        implicit none

        integer, intent(in) :: ord, lb
        real, intent(in) :: dx
        real, dimension(3), intent(in) :: W_0, W_1, W_2, W_3, W_4, W_5, W_6
        real, dimension(3), intent(in) :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12
        real, dimension(3), intent(in) :: E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8
        real, dimension(3), intent(out):: W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7
        real, dimension(6,3) :: tmp_pos
        real, dimension(3) :: xdUM, xdUP, ydUM, ydUP, zdUM, zdUP  
        real, dimension(3) ::  limiter1, limiter, poly, y, sig, dU2, minpath, maxpath
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

        do i= 1, 3
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

        end subroutine field_reconstruct_edge

!-------------------------------------------------------------------------------!   
! Subroutine field_reconstruct_edge1                                            !
!-------------------------------------------------------------------------------!

        subroutine field_reconstruct_edge1(lb, W_0, W_1, W_2, W_3, W_4, W_5, W_6, &
        C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12, E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, &
        W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7, ord,dx)

        implicit none
        
        integer, intent(in) :: ord, lb
        real, intent(in) :: dx
        real, dimension(3), intent(in) :: W_0, W_1, W_2, W_3, W_4, W_5, W_6
        real, dimension(3), intent(in) :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12
        real, dimension(3), intent(in) :: E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8
        real, dimension(3), intent(out):: W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7
        real, dimension(3) :: tmp_pos1, tmp_pos2, tmp_pos3, tmp_pos4, tmp_pos5, tmp_pos6, tmp_pos7, tmp_pos8
        real, dimension(3) :: grad_x, grad_y, grad_z
        real, dimension(3) ::  limiter, minpath, maxpath, mincell, maxcell
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
        
        do i=1,3
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

        end subroutine field_reconstruct_edge1


!-------------------------------------------------------------------------------!   
! Subroutine field_reconstruct_edge2                                            !
!-------------------------------------------------------------------------------!

        subroutine field_reconstruct_edge2(lb, ic, jc, kc, W_0, W_1, W_2, W_3, W_4, W_5, W_6, &
        C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12, E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, &
        W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7, ord)


        implicit none

        integer, intent(in) :: ord, lb, ic, jc, kc
        real, dimension(3), intent(in) :: W_0, W_1, W_2, W_3, W_4, W_5, W_6
        real, dimension(3), intent(in) :: C_1, C_2, C_3, C_4, C_5, C_6, C_7, C_8, C_9, C_10, C_11, C_12
        real, dimension(3), intent(in) :: E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8
        real, dimension(3), intent(out):: W_0_edge1, W_0_edge2, W_0_edge3, W_0_edge4, W_0_edge5, W_0_edge6, W_0_edge7
        real, dimension(26) :: dx, dy, dz
        real, dimension(3,-1:1,-1:1,-1:1) :: Wdat
        real, dimension(3,26) :: dWdat
        real, dimension(3) ::  sum_Wx, sum_Wy, sum_Wz, limiter, minpath, maxpath, mincell, maxcell
        real :: sum_x2, sum_y2, sum_z2, sum_xy, sum_xz, sum_yz, Den
        real, dimension(3) :: tmp_pos1, tmp_pos2, tmp_pos3, tmp_pos4, tmp_pos5, tmp_pos6, tmp_pos7, tmp_pos8
        real, dimension(3) :: grad_x, grad_y, grad_z
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
        
        do i=1,3
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

        end subroutine field_reconstruct_edge2

 end module field_module
