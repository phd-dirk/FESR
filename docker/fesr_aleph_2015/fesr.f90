#/usr/lib/x86_64-linux-gnu
!     Main program for finite energy sum rules (FESR)
!     Last change: Matthias Jamin,  9.10.2014

!-----------------------------------------------------------------------
      program fesr
        use params
        use cVA0_const
        use cVA0_depen

        implicit none
        integer :: npar, iflag = 0
        real(dp) :: chi2, dump, xval(13)
        real(dp), dimension(13) :: grad = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
        complex(dc) :: cintVA0FO, DV0, zarg, DVA4GG, DVA4m4, DVA4qq, cintVA4FO, DVpA68, cintVpA68,cintVpA4FO, deltaP
        real(dp) :: cintDVP0_VA, cintDVP1_VA, cintDVP2_VA, cintDVP3_VA, cintDVP4_VA, cintDV_VA


!     Write out experimental moments of ALEPH spectral functions
!     call writeVAmomAleph()

!     Perform fit of R_tau,V moments
!     call rtauVfit()

!     Perform fit of R_tau,A moments
!     call rtauAfit()

!     Perform fit of R_tau,V and R_tau,A moments simultaneously
!     call rtauVAfit()

!     Perform fit of R_tau,V+A moments
      call rtauVpAfit()

!     Perform fit of R_tau,V-A moments
!     call rtauVmAfit()

!        xval(1) = 0.3179d0 ! astau
!        xval(2) = 2.1d-2 ! aGGInv
!        xval(3) = -0.15d0 ! d6
!        xval(4) = 0.24d0 ! d8
!
!        xval(5) = 2.83d2 ! c_51
!
!        xval(6) = 3.56d0 ! delV
!        xval(7) = 0.58d0 ! gamV
!        xval(8) = -1.92d0 ! alpV
!        xval(9) = 4.07d0 ! betV
!
!        xval(10) = 1.68d0 ! delA
!        xval(11) = 1.41d0 ! gamA
!        xval(12) = 5.16d0 ! alpA
!        xval(13) = 2.13d0 ! betA
!
!        open (unit=7, file='out/minuit_all.out')
!        open (unit=8, file='out/fit_all.out')
!        call chisqVpA(npar, grad, chi2, xval, iflag)
!        print*, "chi2", chi2

!        call init_cVA0nk(5,cV0)
      !function zarg(nf,q2,p2,ap)
!       print*, "amu", zarg(3, (2d0, 0)**2, 3.1572314596d0**2, 0.31927d0/pi)
!        print*, "cV0", cV0
!        print*, "DV0", DV0(3, cV0, 3d0, 3d0, 0.3156d0, 5)


                   ! &( DVA4GG(scon,-smu2,astau,nmax,aGGinv) +&
                   !   &DVA4qq(r,i,j,scon,-smu2,astau,nmax) +&
                   !   &DVA4m4(r,i,j,scon,-smu2,astau,nmax) )
!         print*, "D4", DVA4GG(3d0,3d0,0.3156d0,5,2.1d-2) +&
!                      &DVA4qq(1,1,2,3d0,3d0,0.3156d0,5) +&
!                     &DVA4m4(1,1,2,3d0,3d0,0.3156d0,5)


!       print*, "cint0FO", cintVA0FO(1, cV0, 2d0, 2d0, 0.3156d0, 5)
!       print*, "cintD4FO", cintVA4FO(1, 1,1,2,3d0,3d0,0.3156d0,2.1d-2)
!       print*, "cintD4VpAFO", cintVpA4FO(1,1,2,2d0,2d0,0.3179d0,2.1d-2)

!       print*, "D6FO", DVpA68(1,2,3d0,0.3156d0,-0.5d0,0.1d0)
!        print*, cintVpA68(1,1,2,3d0,0.3156d0,-0.5d0,0.1d0)
         !momth_vplusa(Nsum,Ns0s,s0s,scheme,cV0,ksi,astau,&
         !                    &aGGinv,rhoD6,c8D8,kaV,gaV,alV,beV,&
         !                    &kaA,gaA,alA,beA,momth)

!         print*, deltaP(1, 3d0)
!         print*, deltaP(1, 3d0)
!          print*, "DVP0", cintDVp0_VA(exp(-1.1d0), 1.2d0, 1.3d0, 1.4d0, 3d0)
!          print*, "DVP1", cintDVp1_VA(exp(-1.1d0), 1.2d0, 1.3d0, 1.4d0, 3d0)
!          print*, "DVP2", cintDVp2_VA(exp(-1.1d0), 1.2d0, 1.3d0, 1.4d0, 3d0)
!          print*, "DVP3", cintDVp3_VA(exp(-1.1d0), 1.2d0, 1.3d0, 1.4d0, 3d0)
!          print*, "DVcintVA", cintDV_VA(1,exp(-1.1d0), 1.2d0, 1.3d0, 1.4d0, 3d0)

      stop
      end
!-----------------------------------------------------------------------
