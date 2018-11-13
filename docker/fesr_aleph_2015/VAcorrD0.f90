
!     Routines to provide the vector and axialvector correlators
!     Assumes flavour non-diagonal correlators
!     Additional disconnected fermion loops in case of
!     flavour singlet correlators not yet implemented!!!
!     Last change: Matthias Jamin, 24.1.2010

!-----------------------------------------------------------------------

!     Perturbative Adler function

      function DV0(nf,cV0,s,mu2,astau,nmax)
      use params

      implicit none
      integer,     intent(in) :: nf, nmax
      real(dp),    intent(in) :: cV0(0:12,1:12,3:6), astau
      complex(dc), intent(in) :: s, mu2
      complex(dc)             :: DV0

      integer     :: n, k
      real(dp)    :: atau
      complex(dc) :: L, amu, zarg, sun

      L = log(-s/mu2)
      atau = astau/pi
      amu = zarg(nf,mu2,stau,atau)

      sun = i0

      do n = 1, nmax
         do k = 1, n
            sun = sun + amu**n*k*cV0(n,k,nf)*L**(k-1)
         end do
      end do
      
      DV0 = 1.d0/(4.d0*pi**2)*( cV0(0,1,nf) + sun)

      return
      end function DV0
!-----------------------------------------------------------------------

!     Borel resummed Adler function in Renormalon Pole Model
!     up to factor N_c/(12pi^2) and starting at order alpha_s

      function DV0BS(c51,s,astau)
      use params
      use RPM_const

      implicit none
      real(dp),    intent(in) :: c51
      complex(dc), intent(in) :: s, astau
      complex(dc)             :: DV0BS

      real(dp)    :: b0
      complex(dc) :: atau, as, DV0hat, zarg, RbUV3P, RbIR3P

      atau = astau/pi

      b0 = 9.d0/(4.d0*pi)
      as = zarg(3,s,stau,atau)

      DV0hat = cUV(1,c51)*RbUV3P(2,1,as) +&
              &cIR(2,c51)*RbIR3P(1,2,as) +&
              &cIR(3,c51)*conjg(RbIR3P(1,3,as)) +&
              &cPO(1,c51)*pi*as + cPO(2,c51)*pi**2*b0*as**2

      DV0BS = DV0hat

      return
      end function DV0BS

!-----------------------------------------------------------------------

!     Perturbative vector correlation function Pi_V(s)
!     up to factor N_c/(12pi^2) and starting at order alpha_s

      function PiVhat0(nf,cV0,s,mu2,astau,nmax)
      use params

      implicit none
      integer,     intent(in) :: nf, nmax
      real(dp),    intent(in) :: cV0(0:12,1:12,3:6), astau
      complex(dc), intent(in) :: s, mu2
      complex(dc)             :: PiVhat0

      integer     :: n, k
      real(dp)    :: atau
      complex(dc) :: L, amu, zarg, sun

      L = log(-s/mu2)
      atau = astau/pi
      amu = zarg(nf,mu2,stau,atau)

      sun = i0
      do n = 1, nmax
         do k = 1, n
            sun = sun - amu**n*cV0(n,k,nf)*L**k
         end do
      end do

      PiVhat0 = sun

      return
      end function PiVhat0
!-----------------------------------------------------------------------

!     Perturbative vector spectral function
!     up to factor N_c/(12pi^2) and starting at order alpha_s

      function rhoVhat0(nf,s,mu2,astau,nmax)
      use num_const
      use cVA0_const
      use cVA0_depen

      implicit none
      integer,     intent(in) :: nf, nmax
      real(dp),    intent(in) :: s, mu2, astau
      real(dp)                :: rhoVhat0

      complex(dc) :: PiVhat0

      call init_cVA0nk(nmax,cV0)
      rhoVhat0 = aimag(PiVhat0(nf,cV0,dcmplx(s,0.d0),&
                              &dcmplx(mu2,0.d0),astau,nmax))/pi

      return
      end function rhoVhat0
!-----------------------------------------------------------------------
