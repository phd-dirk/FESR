
!     Compute theoretical moments
!     Last change: Matthias Jamin, 25.6.2011

!-----------------------------------------------------------------------

!     Theoretical moments of the V/A correlators

      subroutine momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,&
                          aGGinv,rhoD6,c8D8,kap,gam,alp,bet,momth)
      use params
      use weights_used

      implicit none
      integer,  intent(in) :: Nsum ! Total number of s_0's
      integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
      real(dp), dimension(Nmom,0:80) :: s0s ! Sets of used s_0's

      character(len=2), intent(in) :: scheme
      integer,  intent(in)  :: r ! r=1 : vector; r=-1 : axialvector
      real(dp), dimension(0:12,1:12,3:6), intent(in) :: cV0
      real(dp), intent(in) :: ksi, astau, aGGinv, rhoD6, c8D8
      real(dp), intent(in) :: kap, gam, alp, bet

      real(dp), intent(out) :: momth(Nsum)

      integer  :: m, n, mn
      real(dp) :: s0, mu2, deltaP
      real(dp) :: cintDVpm_VA !, cintDV_VA
      complex(dc), dimension(Nmom) :: crtauth
      complex(dc) :: cintVA0FO, cintVA0CI, cintVA0BS
      complex(dc) :: cintVA2FO, cintVA4FO, cintVA68
!     complex(dc) :: cintVA2CI, cintVA4CI, cintVA68

!     Compose vector of theoretical moments
      mn = 0
      do m = 1, Nmom
         do n = 1, Ns0s(m)
            mn = mn + 1
            s0  = s0s(m,n)
            mu2 = ksi*s0
            select case (scheme)
            case ("FO")
               crtauth(m) = cintVA0FO(m,cV0,s0,mu2,astau,5)
            case ("CI")
               crtauth(m) = cintVA0CI(m,cV0,s0,mu2,astau,5)
            case ("BS")
               crtauth(m) = cintVA0FO(m,cV0,s0,mu2,astau,0) +&
                           &cintVA0BS(m,cV0(5,1,3),s0,astau)
            end select 
            crtauth(m) = crtauth(m)&
!                     &+ cintVA2FO(m,r,1,2,s0,s0,astau,2)&
!                     &+ cintVA4FO(m,r,1,2,s0,s0,astau,aGGinv)&
                      &+ cintVA68(m,r,1,2,s0,astau,rhoD6,c8D8) !&
!                     &+ cintVA0FO(m,cV0,s0,s0,astau,0)*deltaEW&
!                     &+ cintDV_VA(m,kap,gam,alp,bet,s0)
!   For the moment, depending on the weight, this has to be changed by hand!!!
!                     &+ cintDVpm_VA(m,kap,gam,alp,bet,s0)
            if (r == -1) then
               crtauth(m) = crtauth(m) + 3.d0*dcmplx(deltaP(m,s0),0)
            end if
            momth(mn) = Vud**2*SEW*real(crtauth(m))
         end do
      end do

      return
      end subroutine momth_va
!-----------------------------------------------------------------------

!     Theoretical moments of the V+A correlator

      subroutine momth_vplusa(Nsum,Ns0s,s0s,scheme,cV0,ksi,astau,&
                             &aGGinv,rhoD6,c8D8,kaV,gaV,alV,beV,&
                             &kaA,gaA,alA,beA,momth)
      use params
      use weights_used

      implicit none
      integer,  intent(in) :: Nsum ! Total number of s_0's
      integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
      real(dp), dimension(Nmom,0:80) :: s0s ! Sets of used s_0's

      character(len=2), intent(in) :: scheme
      real(dp), dimension(0:12,1:12,3:6), intent(in) :: cV0
      real(dp), intent(in) :: ksi, astau, aGGinv, rhoD6, c8D8
      real(dp), intent(in) :: kaV, gaV, alV, beV, kaA, gaA, alA, beA

      real(dp), intent(out) :: momth(Nsum)

      integer  :: m, n, mn
      real(dp) :: s0, mu2, deltaP
      real(dp) :: cintDVpm_VA !, cintDV_VA
      complex(dc), dimension(Nmom) :: crtauth
      complex(dc) :: cintVpA0FO, cintVpA0CI, cintVpA0BS
      complex(dc) :: cintVpA2FO, cintVpA4FO, cintVpA68
!     complex(dc) :: cintVpA2CI, cintVpA4CI, cintVpA68

!     Compose vector of theoretical moments
      mn = 0
      do m = 1, Nmom
         do n = 1, Ns0s(m)
            mn = mn + 1
            s0  = s0s(m,n)
            mu2 = ksi*s0
            select case (scheme)
            case ("FO")
               crtauth(m) = cintVpA0FO(m,cV0,s0,mu2,astau,5)
            case ("CI")
               crtauth(m) = cintVpA0CI(m,cV0,s0,mu2,astau,5)
            case ("BS")
               crtauth(m) = cintVpA0FO(m,cV0,s0,mu2,astau,0) +&
                           &cintVpA0BS(m,cV0(5,1,3),s0,astau)
            end select 
            crtauth(m) = crtauth(m)&
!                     &+ cintVpA2FO(m,1,2,s0,s0,astau,2)&
                      &+ cintVpA4FO(m,1,2,s0,s0,astau,aGGinv)&
                      &+ cintVpA68(m,1,2,s0,astau,rhoD6,c8D8)&
!                     &+ cintVpA0FO(m,cV0,s0,s0,astau,0)*deltaEW&
!   For the moment, depending on the weight, this has to be changed by hand!!!
!                     &+ cintDVpm_VA(m,kaV,gaV,alV,beV,s0)&
!                     &+ cintDVpm_VA(m,kaA,gaA,alA,beA,s0)&
                      &+ 3.d0*dcmplx(deltaP(m,s0),0)
            momth(mn) = Vud**2*SEW*real(crtauth(m))
         end do
      end do

      return
      end subroutine momth_vplusa
!-----------------------------------------------------------------------
