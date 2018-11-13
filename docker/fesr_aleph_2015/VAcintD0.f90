
!     Contour integration of the purly perturbative V/A correlators.
!     Results are given for R_tau,V/A^(0) as well as R_tau,V+A^(0),
!     apart from factor |V_ud^2| S_EW, in FOPT, CIPT and for the BS.
!     Last change: Matthias Jamin, 22.1.2010

!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=0) in FOPT

      function cintVpA0FO(m,cV0,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc)          :: cintVpA0FO

      complex(dc)          :: cintVA0FO

      cintVpA0FO = 2.d0*cintVA0FO(m,cV0,s0,mu2,astau,nmax)

      return
      end function cintVpA0FO
!-----------------------------------------------------------------------

!     Contour integral of D_V/A(D=0) in FOPT

      function cintVA0FO(m,cV0,s0,mu2,astau,nmax)

      use params
      use gauss
      use cVA0_depen
      use weights_used

      implicit none
      integer, parameter :: nf = 3

      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc) :: cintVA0FO

      integer     :: n
      complex(dc) :: DV0, scon, sun, x, zmu2

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      sun = i0
      zmu2 = dcmplx(mu2,0)
!     Initialise dependent perturbative coefficients
      call init_cVA0nk(nmax,cV0)

!     R_tau in FOPT
      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x
         sun = sun + wgtD(m,x)*DV0(nf,cV0,scon,zmu2,astau,nmax)*wg(n)
      end do
      
!     R_tau,V/A^(0) in FOPT
      cintVA0FO = 3.d0*pi*sun

      return
      end function cintVA0FO
!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=0) in CIPT

      function cintVpA0CI(m,cV0,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc)          :: cintVpA0CI

      complex(dc)          :: cintVA0CI

      cintVpA0CI = 2.d0*cintVA0CI(m,cV0,s0,mu2,astau,nmax)

      return
      end function cintVpA0CI
!-----------------------------------------------------------------------

!     Contour integral for D_V/A(D=0) in CIPT

      function cintVA0CI(m,cV0,s0,mu2,astau,nmax)
      use params
      use gauss
      use cVA0_depen
      use weights_used

      implicit none
      integer, parameter :: nf = 3

      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc) :: cintVA0CI

      integer     :: n
      complex(dc) :: DV0, scon, smu2, sun, x, zmu2

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      sun = i0
      zmu2 = dcmplx(mu2,0)
!     Initialise dependent perturbative coefficients
      call init_cVA0nk(nmax,cV0)

!     R_tau in CIPT
      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x;  smu2 = zmu2*x
         sun = sun + wgtD(m,x)*DV0(nf,cV0,scon,-smu2,astau,nmax)*wg(n)
      end do
      
!     delta_V+A^(0) in CIPT
      cintVA0CI = 3.d0*pi*sun

      return
      end function cintVA0CI
!-----------------------------------------------------------------------

!     Contour integral for D_V+A(D=0) in BS
!     Without leading-order contribution

      function cintVpA0BS(m,c51,s0,astau)
      use num_const

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: c51, s0, astau
      complex(dc)          :: cintVpA0BS

      complex(dc)          :: cintVA0BS

      cintVpA0BS = 2.d0*cintVA0BS(m,c51,s0,astau)

      return
      end function cintVpA0BS
!-----------------------------------------------------------------------

!     Contour integral for renormalon pole model (RPM)
!     Without leading-order contribution
!!!   c_51 dependence to be implemented !!!

      function cintVA0BS(m,c51,s0,astau)
      use params
      use gauss
      use rge_const
      use weights_used
      use RPM_const

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: c51, s0, astau
      complex(dc)          :: cintVA0BS

      integer     :: n
      complex(dc) :: sumUV1, sumIR2, sumIR3, sumPO1, sumPO2
      complex(dc) :: atau, x, scon, as, zarg
      complex(dc) :: RbUV3P, RbIR3P
      complex(dc) :: conintUV1, conintIR2, conintIR3
      complex(dc) :: conintPO1, conintPO2

!     Initialise Gauss integration
!     call gauleg(-pi,pi,yy,wg,ngau)

      atau = astau/pi

      sumUV1 = i0;
      sumIR2 = i0;  sumIR3 = i0;
      sumPO1 = i0;  sumPO2 = i0;

!     Contour integrals for IR and UV poles as well as polynomial

      do n = 1,ngau
         x = -exp(ii*yy(n));  scon = s0*x
         as = zarg(3,-scon,stau,atau)
         sumUV1 = sumUV1 + wgtD(m,x)*RbUV3P(2,1,as)*wg(n)
         sumIR2 = sumIR2 + wgtD(m,x)*RbIR3P(1,2,as)*wg(n)
         sumIR3 = sumIR3 + wgtD(m,x)*RbIR3P(1,3,as)*wg(n)
         sumPO1 = sumPO1 + wgtD(m,x)*as*wg(n)
         sumPO2 = sumPO2 + wgtD(m,x)*as**2*wg(n)
      end do

      conintUV1 = cUV(1,c51)*sumUV1/(2.d0*pi)

      conintIR2 = cIR(2,c51)*sumIR2/(2.d0*pi)
      conintIR2 = cmplx(real(conintIR2),abs(aimag(conintIR2)),dc)
      conintIR3 = cIR(3,c51)*sumIR3/(2.d0*pi)
      conintIR3 = cmplx(real(conintIR3),abs(aimag(conintIR3)),dc)

      conintPO1 = cPO(1,c51)*pi*sumPO1/(2.d0*pi)
      conintPO2 = cPO(2,c51)*pi*(be_1(3)/2.d0)*sumPO2/(2.d0*pi)

!     write (*,*) 'conintIR2 = ', conintIR2
!     write (*,*) 'conintIR3 = ', conintIR3
!     write (*,*) 'conintUV1 = ', conintUV1
!     write (*,*) 'conintPO  = ', conintPO1 + conintPO2

!     delta_V+A^(0) in RPM
      cintVA0BS = 1.5d0*( conintUV1 + conintIR2 + conintIR3 +&
                         &conintPO1 + conintPO2 )

      return
      end function cintVA0BS
!-----------------------------------------------------------------------
