
!     Contour integration of the D=6,8 contribution
!     Last change: Matthias Jamin, 22.1.2010

!-----------------------------------------------------------------------

!     Contour integral of D_V/A(D=6+8)

      function cintVA68(m,r,i,j,s0,astau,rhoVA,c8VA)
      use params
      use gauss
      use weights_used

      implicit none
      integer,  intent(in) :: m, r, i, j
      real(dp), intent(in) :: s0, astau, rhoVA, c8VA
      complex(dc)          :: cintVA68

      integer     :: n
      complex(dc) :: scon, sun, x, DVA68

      sun = i0

      do n = 1, ngau
         x = -exp(ii*yy(n));  scon = s0*x
         sun = sun + wgtD(m,x)*DVA68(r,i,j,scon,astau,rhoVA,c8VA)*wg(n)
      end do

!     Corresponds to delta_V/A^(6+8)
      cintVA68 = 3.d0*pi*sun

      return
      end function cintVA68
!-----------------------------------------------------------------------

!     Contour integral of D_V+A(D=6+8)

      function cintVpA68(m,i,j,s0,astau,rhoVpA,c8VpA)
      use params
      use gauss
      use weights_used

      implicit none
      integer,  intent(in) :: m, i, j
      real(dp), intent(in) :: s0, astau, rhoVpA, c8VpA
      complex(dc)          :: cintVpA68

      integer     :: n
      complex(dc) :: scon, sun, x, DVpA68

      sun = i0

      do n = 1, ngau
         x = -exp(ii*yy(n));  scon = s0*x
         sun = sun + wgtD(m,x)*DVpA68(i,j,scon,astau,rhoVpA,c8VpA)*wg(n)
      end do

!     Corresponds to delta_V+A^(6+8)
      cintVpA68 = 3.d0*pi*sun

      return
      end function cintVpA68
!-----------------------------------------------------------------------
