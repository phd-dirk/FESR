
!     Contribution of the pion pole to axialvector, V+A and V-A
!     Last change: Matthias Jamin, 19.3.2011

!-----------------------------------------------------------------------

!     R_tau,A(s_0, pion pole)^w

      subroutine pimom(Nsum,Ns0s,s0s,mom,pifac,dpifac)
      use params
      use weights_used

      implicit none
      integer,  intent(in) :: Nsum ! Total number of s_0's
      integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
      real(dp), dimension(Nmom,0:80) :: s0s ! Sets of used s_0's
      real(dp), intent(out) :: mom(Nsum), pifac, dpifac

      integer  :: k, l, kl
      real(dp) :: momA(Nsum), momP(Nsum)

       pifac = 24.d0*(pi*Vud*fpi)**2*SEW
      dpifac = pifac*sqrt(4.d0*(dVud/Vud)**2+(dSEW/SEW)**2+&
                         &4.d0*(dfpi/fpi)**2)

      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
!     Axialvector pion-pole contribution
            momA(kl) = pifac/s0s(k,l)*&
                      &real(wgtr(k,dcmplx(mpim**2/s0s(k,l),0)))
!     Pseudoscalar pion-pole contribution
            momP(kl) = momA(kl)*(-2.d0*mpim**2/(stau+2.d0*mpim**2))
            mom(kl) = momA(kl) + momP(kl)
         end do
      end do

      return
      end subroutine pimom
!-----------------------------------------------------------------------
