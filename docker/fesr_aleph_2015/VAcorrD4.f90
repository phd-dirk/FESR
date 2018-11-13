
!     Routine to provide the vector and axialvector correlators at D=4
!     Only for THREE quark flavours!
!     Last change: Matthias Jamin, 28.12.2009

!-----------------------------------------------------------------------

!     Gluon condensate contribution

      function DVA4GG(s,mu2,astau,nmax,aGGinv)
      use params

      implicit none

      integer,     intent(in) :: nmax
      real(dp),    intent(in) :: astau, aGGinv
      complex(dc), intent(in) :: s, mu2
      complex(dc)             :: DVA4GG

      real(dp)    :: atau
!     This coefficient is not yet know, but the result should
!     not depend on it. So we set it to zero.
      real(dp)    :: pLT3 = 0.d0
      complex(dc) :: amu, zarg, sun

      atau = astau/pi
      amu = zarg(3,mu2,stau,atau)

      sun = i0

!     Sum consecutive orders depending on n_max

      if (nmax>0) then
         sun = 1.d0/6.d0
      end if
      if (nmax>1) then
         sun = sun - 11.d0/108.d0*amu
      end if
      if (nmax>2) then
         sun = sun + (11.d0/48.d0*log(-s/mu2)+pLT3/6.d0-&
                     &8773.d0/15552.d0)*amu**2
      end if

      DVA4GG = sun*aGGinv/s**2

      return
      end function DVA4GG
!-----------------------------------------------------------------------

!     Quark condensate contribution

      function DVA4qq(r,i,j,s,mu2,astau,nmax)
      use condens

      implicit none

      integer,     intent(in) :: r, i, j, nmax
      real(dp),    intent(in) :: astau
      complex(dc), intent(in) :: s, mu2
      complex(dc)             :: DVA4qq

      real(dp)    :: atau, mqqa, mqqb, mqqs
!     These coefficients are not yet know, but the result should
!     not depend on them. So we set them to zero.
      real(dp)    :: pLT3 = 0.d0, qLT3 = 0.d0, rLT3 = 0.d0, tLT3 = 0.d0
      complex(dc) :: L, amu, zarg, sun

      L = log(-s/mu2)
      atau = astau/pi
      amu = zarg(3,mu2,stau,atau)

      mqqa = mq(i)*qqinv(i) + mq(j)*qqinv(j)
      mqqb = r*(mq(i)*qqinv(j) + mq(j)*qqinv(i))
      mqqs = mq(1)*qqinv(1) + mq(2)*qqinv(2) + mq(3)*qqinv(3)

      sun = i0

!     Sum consecutive orders depending on n_max

      if (nmax>-1) then
         sun = 2.d0*mqqa
      end if
      if (nmax>0) then
         sun = sun + (-2.d0*mqqa+8.d0/3.d0*mqqb+8.d0/27.d0*mqqs)*amu
      end if
      if (nmax>1) then
         sun = sun + ((4.5d0*L-131.d0/12.d0)*mqqa+(-6.d0*L+68.d0/3.d0)*&
                     &mqqb+(-2.d0/3.d0*L-176.d0/243.d0+8.d0/3.d0*z3)*&
                     &mqqs)*amu**2
      end if
      if (nmax>2) then
         sun = sun + ((-81.d0/8.d0*L*L+457.d0/8.d0*L+2.d0*qLT3)*mqqa+&
                     &(27.d0/2.d0*L*L-338.d0/3.d0*L+8.d0/3.d0*tLT3)*&
                     &mqqb+(1.5d0*L*L+(56.d0/27.d0-12.d0*z3)*L+50407/&
                           &17496d0+8.d0/27.d0*pLT3+rLT3-20.d0/27.d0*z3&
                           &)*mqqs)*amu**3
      end if

      DVA4qq = sun/s**2

      return
      end function DVA4qq
!-----------------------------------------------------------------------

!     m_q^4 contribution
!     Constant terms at order a[mu^2] not implemented

      function DVA4m4(r,i,j,s,mu2,astau,nmax)
      use condens

      implicit none

      integer,     intent(in) :: r, i, j, nmax
      real(dp),    intent(in) :: astau
      complex(dc), intent(in) :: s, mu2
      complex(dc)             :: DVA4m4

      real(dp)    :: atau
      complex(dc) :: L, amu, rmq, zarg, zmrg, m4a, m4b, m4c, m4d, sun

      L = log(-s/mu2)
      atau = astau/pi
      amu = zarg(3,mu2,stau,atau)
      rmq = zmrg(3,mu2,stau,atau)

!     Different mass combinations at scale mu^2
      m4a = (mq(i)**4 + mq(j)**4)*rmq**4
      m4b = r*(mq(i)*mq(j)**3 + mq(j)*mq(i)**3)*rmq**4
      m4c = (mq(i)**2*mq(j)**2)*rmq**4
      m4d = (mq(1)**4 + mq(2)**4 + mq(3)**4)*rmq**4

      sun = i0

!     Sum consecutive orders depending on n_max

      if (nmax>1) then
         sun = sun + ( - 6.d0/7.d0/amu + 1.5d0*L - 0.25d0 )*m4a -&
                    &8.d0/7.d0*m4b + 3.d0*m4c - m4d/14.d0

      end if
      if (nmax>2) then
         sun = sun + ( - 3.d0*L*L + 74.d0/7.d0*L )*amu*m4a +&
                    &32.d0/7.d0*L*amu*m4b - 12.d0*L*amu*m4c +&
                    &2.d0/7.d0*L*amu*m4d
      end if

      DVA4m4 = sun/(pi*s)**2

      return
      end function DVA4m4
!-----------------------------------------------------------------------
