
!     Contour integration of the D=2 contribution
!     Last change: Matthias Jamin, 22.1.2010

!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=2) in FOPT

      function cintVpA2FO(m,i,j,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: cintVpA2FO

      complex(dc)          :: cintVA2FO

      cintVpA2FO = cintVA2FO(m, 1,i,j,s0,mu2,astau,nmax) +&
                  &cintVA2FO(m,-1,i,j,s0,mu2,astau,nmax)

      return
      end function cintVpA2FO
!-----------------------------------------------------------------------

!     Contour integral of D_V/A(D=2) in FOPT

      function cintVA2FO(m,r,i,j,s0,mu2,astau,nmax)
      use params
      use gauss
      use weights_used

      implicit none
      integer,  intent(in) :: m, r, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: cintVA2FO

      integer     :: n
      complex(dc) :: scon, sun, x, zmu2, DVA2m2

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      sun = i0
      zmu2 = dcmplx(mu2,0)

      sun = i0

      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x
         sun = sun + wgtD(m,x)*DVA2m2(r,i,j,scon,zmu2,astau,nmax)*wg(n)
      end do

!     Corresponds to delta_V/A^(2)
      cintVA2FO = 3.d0*pi*sun

      return
      end function cintVA2FO
!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=2) in CIPT

      function cintVpA2CI(m,i,j,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: cintVpA2CI

      complex(dc)          :: cintVA2CI

      cintVpA2CI = cintVA2CI(m, 1,i,j,s0,mu2,astau,nmax) +&
                  &cintVA2CI(m,-1,i,j,s0,mu2,astau,nmax)

      return
      end function cintVpA2CI
!-----------------------------------------------------------------------

!     Contour integral of D_V/A(D=2) in CIPT

      function cintVA2CI(m,r,i,j,s0,mu2,astau,nmax)
      use params
      use gauss
      use weights_used

      implicit none
      integer,  intent(in) :: m, r, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: cintVA2CI

      integer     :: n
      complex(dc) :: scon, smu2, sun, x, zmu2, DVA2m2

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      sun = i0
      zmu2 = dcmplx(mu2,0)

      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x;  smu2 = zmu2*x
         sun = sun + wgtD(m,x)*DVA2m2(r,i,j,scon,-smu2,astau,nmax)*wg(n)
      end do

!     Corresponds to delta_V/A^(2)
      cintVA2CI = 3.d0*pi*sun

      return
      end function cintVA2CI
!-----------------------------------------------------------------------
