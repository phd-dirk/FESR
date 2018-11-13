
!     Routine to provide the vector and axialvector correlators at D=6,8
!     Only for THREE quark flavours!
!     Last change: Matthias Jamin, 21.9.2011

!-----------------------------------------------------------------------

!     V/A correlator at D=6,8

      function DVA68(r,i,j,s,astau,rhoVA,c8VA)
      use condens

      implicit none
      integer,  intent(in) :: r, i, j
      real(dp),    intent(in) :: astau, rhoVA, c8VA
      complex(dc), intent(in) :: s
      complex(dc)             :: DVA68

!     OLD version
!     DVA68 = r*32.d0/3.d0*pi*astau*rhoVA*qqmtau(i)*qqmtau(j)/s**3 -&
!            &32.d0/27.d0*pi*astau*rhoVA*(qqmtau(i)**2 + qqmtau(j)**2)/&
!            &s**3 + c8VA/s**4

!     select case (r)
!     case ( 1)
!        DVA68 = 32.d0/27.d0*pi*astau*7.d0*rhoVA*qqmtau(1)**2/s**3+&
!               &4.d-2*c8VA/s**4
!     case (-1)
!        DVA68 = 32.d0/27.d0*pi*astau*(-11.d0)*rhoVA*qqmtau(1)**2/s**3+&
!               &4.d-2*c8VA/s**4
!     end select

!     Parametrisation of article
!     rhoVA = - 10^2 C_6V/A;  c8VA = 10^2 C_8V/A      
      DVA68 = 3.d-2*rhoVA/s**3 + 4.d-2*c8VA/s**4

      return
      end function DVA68
!-----------------------------------------------------------------------

!     V+A correlator at D=6,8

      function DVpA68(i,j,s,astau,rhoVpA,c8VpA)
      use condens

      implicit none
      integer,  intent(in) :: i, j
      real(dp),    intent(in) :: astau, rhoVpA, c8VpA
      complex(dc), intent(in) :: s
      complex(dc)             :: DVpA68

!     DVpA68 = - 64.d0/27.d0*pi*astau*rhoVpA*(qqmtau(i)**2 +&
!               &qqmtau(j)**2)/s**3 + 2.d0*c8VpA/s**4

      DVpA68 = 3.d-2*rhoVpA/s**3 + 4.d-2*c8VpA/s**4

      return
      end function DVpA68
!-----------------------------------------------------------------------
