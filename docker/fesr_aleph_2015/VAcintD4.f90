
!     Contour integration of the D=4 contribution
!     Last change: Matthias Jamin, 22.1.2010

!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=4) in FOPT

      function cintVpA4FO(m,i,j,s0,mu2,astau,aGGinv)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: cintVpA4FO

      complex(dc)          :: cintVA4FO

      cintVpA4FO = cintVA4FO(m, 1,i,j,s0,mu2,astau,aGGinv) +&
                  &cintVA4FO(m,-1,i,j,s0,mu2,astau,aGGinv)

      return
      end function cintVpA4FO
!-----------------------------------------------------------------------

!     Contour integral of D_V/A(D=4) in FOPT

      function cintVA4FO(m,r,i,j,s0,mu2,astau,aGGinv)
      use params
      use gauss
      use weights_used

      implicit none
      integer,  intent(in) :: m, r, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: cintVA4FO

      integer     :: n, nmax
      complex(dc) :: scon, sun, x, zmu2, DVA4GG, DVA4qq, DVA4m4

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      sun = i0
      zmu2 = dcmplx(mu2,0)

!     Test if the pure 1/s^2 term integrates to zero
      do n = 1,ngau
         x = -exp(ii*yy(n))
         sun = sun + wgtD(m,x)/x**2*wg(n)
      end do

      nmax = 2
      if (abs(real(sun))<2.d-14) then
         nmax = 3
      end if

      sun = i0

      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x
         sun = sun + wgtD(m,x)*wg(n)*&
                    &( DVA4GG(scon,zmu2,astau,nmax,aGGinv) +&
                      &DVA4qq(r,i,j,scon,zmu2,astau,nmax) +&
                      &DVA4m4(r,i,j,scon,zmu2,astau,nmax) )
      end do

!     Corresponds to delta_V/A^(4)
      cintVA4FO = 3.d0*pi*sun

      return
      end function cintVA4FO
!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=4) in CIPT

      function cintVpA4CI(m,i,j,s0,mu2,astau,aGGinv)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: cintVpA4CI

      complex(dc)          :: cintVA4CI

      cintVpA4CI = cintVA4CI(m, 1,i,j,s0,mu2,astau,aGGinv) +&
                  &cintVA4CI(m,-1,i,j,s0,mu2,astau,aGGinv)

      return
      end function cintVpA4CI
!-----------------------------------------------------------------------

!     Contour integral in CIPT

      function cintVA4CI(m,r,i,j,s0,mu2,astau,aGGinv)
      use params
      use gauss
      use weights_used

      implicit none
      integer,  intent(in) :: m, r, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: cintVA4CI

      integer     :: n, nmax = 2 ! Only for CIPT
      complex(dc) :: scon, smu2, sun, x, zmu2, DVA4GG, DVA4qq, DVA4m4

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      zmu2 = dcmplx(mu2,0)
      sun = i0

      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x;  smu2 = zmu2*x
         sun = sun + wgtD(m,x)*wg(n)*&
                    &( DVA4GG(scon,-smu2,astau,nmax,aGGinv) +&
                      &DVA4qq(r,i,j,scon,-smu2,astau,nmax) +&
                      &DVA4m4(r,i,j,scon,-smu2,astau,nmax) )
      end do

!     Corresponds to delta_V/A^(4)
      cintVA4CI = 3.d0*pi*sun

      return
      end function cintVA4CI
!-----------------------------------------------------------------------
