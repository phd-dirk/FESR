
!     Duality violating contribution for V/A according
!     to the model by Cata, Golterman and Peris
!     Last change: Matthias Jamin, 24.5.2010

!-----------------------------------------------------------------------

!     1/pi Im rho_DV,V(s)

      function rhoDV_VA(ka,ga,al,be,s)
      use num_const

      implicit none
      real(dp), intent(in) :: ka, ga, al, be, s
      real(dp)             :: rhoDV_VA

      rhoDV_VA = ka*exp(-ga*s)*sin(al+be*s)

      return
      end function rhoDV_VA
!-----------------------------------------------------------------------

!     Moment integral for duality violating contribution

      function cintDV_VA(m,ka,ga,al,be,s0)
      use num_const
      use weights_used

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: ka, ga, al, be, s0
      real(dp)             :: cintDV_VA

      integer :: n
      integer,  parameter :: ngau = 1201
      real(dp), parameter :: bb = -5.d0
      real(dp) :: s, sun, rhoDV_VA
      real(dp), dimension(ngau) :: zz, wz

!     Initialise Gauss integration
      call gauleg(0.d0,1.d0,zz,wz,ngau)

      sun = 0.d0

      do n = 1, ngau
      s = s0*(1.d0-bb*zz(n))/(1.d0-zz(n))
      sun = sun + (1.d0-bb)/(1.d0-zz(n))**2*&
                 &real(wgtr(m,dcmplx(s/s0,0)))*&
                 &rhoDV_VA(ka,ga,al,be,s)*wz(n)
      end do

      cintDV_VA = - 12.d0*pi**2*sun

      return
      end function cintDV_VA
!-----------------------------------------------------------------------

!     Moment integral for duality violating contribution
!     with power-like weight function w(x) = x^m

      function cintDVpm_VA(m,ka,ga,al,be,s0)
      use num_const
      use weights_used

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: ka, ga, al, be, s0
      real(dp)             :: cintDVpm_VA

      real(dp)             :: cintDVp0_VA, cintDVp1_VA,&
                             &cintDVp2_VA, cintDVp3_VA

      cintDVpm_VA = 0.d0

      select case (m)
      case (1)
!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0)
         cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
                 &3.d0*cintDVp2_VA(ka,ga,al,be,s0) +&
                 &2.d0*cintDVp3_VA(ka,ga,al,be,s0)

!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
!                &3.d0*cintDVp1_VA(ka,ga,al,be,s0) +&
!                &3.d0*cintDVp2_VA(ka,ga,al,be,s0) -&
!                     &cintDVp3_VA(ka,ga,al,be,s0)
      case (2)
!        cintDVpm_VA = cintDVp1_VA(ka,ga,al,be,s0)
!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
!                     &cintDVp1_VA(ka,ga,al,be,s0)
!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
!                     &cintDVp2_VA(ka,ga,al,be,s0)
      case (3)
!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
!                &3.d0*cintDVp2_VA(ka,ga,al,be,s0) +&
!                &2.d0*cintDVp3_VA(ka,ga,al,be,s0)
      case (4)
!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
!                     &cintDVp3_VA(ka,ga,al,be,s0)
      case (5)
!        cintDVpm_VA = cintDVp0_VA(ka,ga,al,be,s0) -&
!                &3.d0*cintDVp2_VA(ka,ga,al,be,s0) +&
!                &2.d0*cintDVp3_VA(ka,ga,al,be,s0)
      end select

      return
      end function cintDVpm_VA
!-----------------------------------------------------------------------

!     Moment integral for duality violating contribution
!     with power-like weight function w(x) = 1

      function cintDVp0_VA(ka,ga,al,be,s0)
      use num_const
      use weights_used

      implicit none
      real(dp), intent(in) :: ka, ga, al, be, s0
      real(dp)             :: cintDVp0_VA

      cintDVp0_VA = - 12.d0*pi**2*ka/(s0*(be*be+ga*ga)*exp(ga*s0))*&
                     &(be*cos(al+be*s0) + ga*sin(al+be*s0))

      return
      end function cintDVp0_VA
!-----------------------------------------------------------------------

!     Moment integral for duality violating contribution
!     with power-like weight function w(x) = x

      function cintDVp1_VA(ka,ga,al,be,s0)
      use num_const
      use weights_used

      implicit none
      real(dp), intent(in) :: ka, ga, al, be, s0
      real(dp)             :: cintDVp1_VA

      real(dp) :: be2, ga2

      be2 = be*be;  ga2 = ga*ga

      cintDVp1_VA = - 12.d0*pi**2*ka/(s0*(be2+ga2))**2/exp(ga*s0)*&
                   &(be*(2.d0*ga+(be2+ga2)*s0)*cos(al+be*s0) +&
                    &(ga2*(1.d0+ga*s0)-be2*(1.d0-ga*s0))*sin(al+be*s0))

      return
      end function cintDVp1_VA
!-----------------------------------------------------------------------

!     Moment integral for duality violating contribution
!     with power-like weight function w(x) = x^2

      function cintDVp2_VA(ka,ga,al,be,s0)
      use num_const
      use weights_used

      implicit none
      real(dp), intent(in) :: ka, ga, al, be, s0
      real(dp)             :: cintDVp2_VA

      real(dp) :: be2, ga2, be2ga2

      be2 = be*be;  ga2 = ga*ga;  be2ga2 = be2+ga2

      cintDVp2_VA = - 12.d0*pi**2*ka/(s0*be2ga2)**3/exp(ga*s0)*&
        &(be*(6.d0*ga2-2.d0*be2+4.d0*ga*be2ga2*s0+(be2ga2*s0)**2)*&
        &cos(al+be*s0) + (2.d0*ga*(ga2-3.d0*be2)+2.d0*(ga**4-be**4)*s0+&
        &ga*(be2ga2*s0)**2)*sin(al+be*s0))

      return
      end function cintDVp2_VA
!-----------------------------------------------------------------------

!     Moment integral for duality violating contribution
!     with power-like weight function w(x) = x^3

      function cintDVp3_VA(ka,ga,al,be,s0)
      use num_const
      use weights_used

      implicit none
      real(dp), intent(in) :: ka, ga, al, be, s0
      real(dp)             :: cintDVp3_VA

      real(dp) :: be2, ga2, be2ga2

      be2 = be*be;  ga2 = ga*ga;  be2ga2 = be2+ga2

      cintDVp3_VA = - 12.d0*pi**2*ka/(s0*be2ga2)**4/exp(ga*s0)*&

        &(be*(24.d0*ga*(ga2-be2)-6.d0*(be2-3.d0*ga2)*be2ga2*s0+&
         &6.d0*ga*(be2ga2*s0)**2+(be2ga2*s0)**3)*cos(al+be*s0) +&
          (be**6*s0**2*(ga*s0-3.d0)+3.d0*be**4*(2.d0+ga*s0*&
          &(2.d0+ga*s0)*(ga*s0-3.d0))+3.d0*(be*ga)**2*(ga*s0*& 
          &(ga*s0*(1.d0+ga*s0)-4.d0)-12.d0)+ga**4*(6.d0+ga*s0*&
          &(6.d0+ga*s0*(3.d0+ga*s0))))*sin(al+be*s0))

      return
      end function cintDVp3_VA
!-----------------------------------------------------------------------
