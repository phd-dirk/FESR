
!     Functions for Renormalon Pole Model (RPM)
!     Last change: Matthias Jamin, 5.6.2009

!-----------------------------------------------------------------------

!     Three IR pole model to include 4-loop running

      function RbIR3P(k,p,a)
      use RPM_const

      implicit none
      integer,     intent(in) :: k, p
      complex(dc), intent(in) :: a
      complex(dc)             :: RbIR3P

      real(dp)    :: ga, bPIR_1, bPIR_2
      complex(dc) :: RbIR1P

      ga = gaIR(k,p)

      RbIR3P = RbIR1P(p,ga,a) + bPIR_1(k,p)*RbIR1P(p,ga-1.d0,a)&
                             &+ bPIR_2(k,p)*RbIR1P(p,ga-2.d0,a)

      return
      end function RbIR3P

!-----------------------------------------------------------------------

!     Single IR pole of general structure

      function RbIR1P(p,ga,a)
      use num_const

      implicit none
      integer,     intent(in) :: p
      real(dp),    intent(in) :: ga
      complex(dc), intent(in) :: a
      complex(dc)             :: RbIR1P

      integer     :: s
      real(dp)    :: b0
      complex(dc) :: z, expintE

      s = 1  !  Sign in complex ambiguity

      b0 = 9.d0/(4.d0*pi)
      z  = p/(pi*b0*a)

      RbIR1P = pi*a/(pi*b0*a)**ga*exp(-z)*&
              &( - z**(1.d0-ga)*expintE(ga,-z) +&
                 ( exp(ii*pi*s*ga) - exp(ii*pi*sign(1.d0,aimag(a))*ga)&
                  &)*dgamma(1.d0-ga) )

      return
      end function RbIR1P

!-----------------------------------------------------------------------

      function bPIR_1(k,p)
      use rge_const
      use RPM_const

      implicit none
      integer, intent(in) :: k, p
      real(dp)            :: bPIR_1

      bPIR_1 = 2.d0*(bb_1(p)+cc_1(p))/be_1(3)/(gaIR(k,p)-1.d0)

      return
      end function bPIR_1

!-----------------------------------------------------------------------

      function bPIR_2(k,p)
      use rge_const
      use RPM_const

      implicit none
      integer, intent(in) :: k, p
      real(dp)            :: bPIR_2

      bPIR_2 = 4.d0*(bb_2(p)+bb_1(p)*cc_1(p)+cc_2(p))/be_1(3)**2/&
              &(gaIR(k,p)-1.d0)/(gaIR(k,p)-2.d0)

      return
      end function bPIR_2

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!     Three UV pole model to include 4-loop running

      function RbUV3P(k,p,a)
      use RPM_const

      implicit none
      integer,     intent(in) :: k, p
      complex(dc), intent(in) :: a
      complex(dc)             :: RbUV3P

      real(dp)    :: ga, bPUV_1, bPUV_2
      complex(dc) :: RbUV1P

      ga = gaUV(k,p)

      RbUV3P = RbUV1P(p,ga,a) + bPUV_1(k,p)*RbUV1P(p,ga-1.d0,a)&
                             &+ bPUV_2(k,p)*RbUV1P(p,ga-2.d0,a)

      return
      end function RbUV3P

!-----------------------------------------------------------------------

!     Single UV pole of general structure

      function RbUV1P(p,ga,a)
      use num_const

      implicit none
      integer,     intent(in) :: p
      real(dp),    intent(in) :: ga
      complex(dc), intent(in) :: a
      complex(dc)             :: RbUV1P

      real(dp)    :: b0
      complex(dc) :: z, incgam

      b0 = 9.d0/(4.d0*pi)
      z  = p/(pi*b0*a)

      RbUV1P = pi*a/(pi*b0*a)**ga*exp(z)*incgam(1.d0-ga,z)

      return
      end function RbUV1P

!-----------------------------------------------------------------------

      function bPUV_1(k,p)
      use num_const

      implicit none
      integer, intent(in) :: k, p
      real(dp)            :: bPUV_1

      real(dp) :: bPIR_1

      bPUV_1 = - bPIR_1(k,-p)

      return
      end function bPUV_1

!-----------------------------------------------------------------------

      function bPUV_2(k,p)
      use num_const

      implicit none
      integer, intent(in) :: k, p
      real(dp)            :: bPUV_2

      real(dp) :: bPIR_2

      bPUV_2 = bPIR_2(k,-p)

      return
      end function bPUV_2

!-----------------------------------------------------------------------
