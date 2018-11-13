
!     Compute derivatives of theoretical moments
!     with respect to the theoretical parameters
!     Last change: Matthias Jamin, 24.2.2012

!-----------------------------------------------------------------------

!     Derivatives of theoretical moments of the V/A correlators

      subroutine dmomth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,&
                           aGGinv,rhoD6,c8D8,kap,gam,alp,bet,dmomth)
      use num_const
      use weights_used

      implicit none
      integer,  intent(in) :: Nsum ! Total number of s_0's
      integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
      real(dp), dimension(Nmom,0:140) :: s0s ! Sets of used s_0's

      character(len=2), intent(in) :: scheme
      integer,  intent(in)  :: r ! r=1 : vector; r=-1 : axialvector
      real(dp), dimension(0:12,1:12,3:6), intent(in) :: cV0
      real(dp), intent(in) :: ksi, astau, aGGinv, rhoD6, c8D8
      real(dp), intent(in) :: kap, gam, alp, bet

!     Derivative with respect to aGGinv not implemented
!     Thus, there are 7 theoretical parameters per channel
      real(dp), intent(out) :: dmomth(7,Nsum)

      integer :: n

      real(dp), parameter :: oep = 1.d-4, tep = 2.d-4
      real(dp) :: momthpoep(Nsum), momthmoep(Nsum),&
                 &momthptep(Nsum), momthmtep(Nsum)

!     Derivative with respect to alpha_s
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau+oep,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau-oep,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau+tep,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau-tep,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet,momthmtep)

      do n = 1, Nsum
         dmomth(1,n) = (8.d0*(momthpoep(n)-momthmoep(n))-&
                       &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

!     Derivative with respect to rhoD6
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6+oep,c8D8,kap,gam,alp,bet,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6-oep,c8D8,kap,gam,alp,bet,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6+tep,c8D8,kap,gam,alp,bet,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6-tep,c8D8,kap,gam,alp,bet,momthmtep)

      do n = 1, Nsum
         dmomth(2,n) = (8.d0*(momthpoep(n)-momthmoep(n))-&
                       &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

!     Derivative with respect to c8D8
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8+oep,kap,gam,alp,bet,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8-oep,kap,gam,alp,bet,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8+tep,kap,gam,alp,bet,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8-tep,kap,gam,alp,bet,momthmtep)

      do n = 1, Nsum
         dmomth(3,n) = (8.d0*(momthpoep(n)-momthmoep(n))-&
                       &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

!     Derivative with respect to delta
!     Be aware of the replacement of delta -> kappa
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap+oep,gam,alp,bet,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap-oep,gam,alp,bet,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap+tep,gam,alp,bet,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap-tep,gam,alp,bet,momthmtep)

      do n = 1, Nsum
         dmomth(4,n) = -kap*(8.d0*(momthpoep(n)-momthmoep(n))-&
                            &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

!     Derivative with respect to gamma
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam+oep,alp,bet,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam-oep,alp,bet,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam+tep,alp,bet,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam-tep,alp,bet,momthmtep)

      do n = 1, Nsum
         dmomth(5,n) = (8.d0*(momthpoep(n)-momthmoep(n))-&
                       &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

!     Derivative with respect to alpha
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp+oep,bet,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp-oep,bet,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp+tep,bet,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp-tep,bet,momthmtep)

      do n = 1, Nsum
         dmomth(6,n) = (8.d0*(momthpoep(n)-momthmoep(n))-&
                       &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

!     Derivative with respect to beta
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet+oep,momthpoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet-oep,momthmoep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet+tep,momthptep)
      call momth_va(Nsum,Ns0s,s0s,scheme,r,cV0,ksi,astau,aGGinv,&
                   &rhoD6,c8D8,kap,gam,alp,bet-tep,momthmtep)

      do n = 1, Nsum
         dmomth(7,n) = (8.d0*(momthpoep(n)-momthmoep(n))-&
                       &momthptep(n)+momthmtep(n))/(6.d0*tep)
      end do

      return
      end subroutine dmomth_va
!-----------------------------------------------------------------------
