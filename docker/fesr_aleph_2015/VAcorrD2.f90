
!     Routine to provide the vector and axialvector correlators at D=2
!     Only for THREE quark flavours!
!     Last change: Matthias Jamin, 22.1.2010

!-----------------------------------------------------------------------

!     m^2 correction to V/A correlators

      function DVA2m2(r,i,j,s,mu2,astau,nmax)
      use params

      implicit none

      integer,     intent(in) :: r, i, j, nmax
      real(dp),    intent(in) :: astau
      complex(dc), intent(in) :: s, mu2
      complex(dc)             :: DVA2m2

      real(dp)    :: atau, ePLT3 = 1.d2 ! Guesstimate
      complex(dc) :: L, amu, rmq, zarg, zmrg, m2a, m2b, m2c, sun

      L = log(-s/mu2)
      atau = astau/pi
      amu = zarg(3,mu2,stau,atau)
      rmq = zmrg(3,mu2,stau,atau)

!     Different mass combinations at scale mu^2
      m2a = (mq(i)**2 + mq(j)**2)*rmq**2
      m2b = r*mq(i)*mq(j)*rmq**2
      m2c = (mq(1)**2 + mq(2)**2 + mq(3)**2)*rmq**2

      sun = i0

!     Sum consecutive orders depending on n_max

      if (nmax>-1) then
         sun = sun + m2a
      end if
      if (nmax>0) then
         sun = sun + ((13.d0/3.d0-2.d0*L)*m2a+2.d0/3.d0*m2b)*amu
      end if
      if (nmax>1) then
         sun = sun + ((4.25d0*L*L-26.d0*L+23077.d0/432.d0+&
                      &179.d0/54.d0*z3-520.d0/27.d0*z5)*m2a +&
                     &(-17.d0/6.d0*L+769.d0/54.d0-55.d0/27.d0*z3-&
                      &5.d0/27.d0*z5)*m2b +&
                     &(-32.d0/9.d0+8.d0/3.d0*z3)*m2c)*amu**2
      end if
      if (nmax>2) then
         sun = sun + ((-221.d0/24.d0*L**3+1153.d0/12.d0*L*L+&
      &(-46253.d0/108.d0-1787.d0/108.d0*z3+3380.d0/27.d0*z5)*L+&
      &3909929.d0/5184.d0-pi**4/36.d0-1541.d0/648.d0*z3+26.5d0*z3**2-&
      &54265.d0/108.d0*z5+79835.d0/648.d0*z7)*m2a +&
                     &(221.d0/24.d0*L*L+(-10831.d0/108.d0+715.d0/54.d0*&
      &z3+65.d0/54.d0*z5)*L+4421.d0/54+ePLT3-715.d0/54.d0*z3-&
      &65.d0/54.d0*z5)*m2b +&
                     &((208.d0/9.d0-52.d0/3.d0*z3)*L-2222.d0/27.d0+&
      &1592.d0/27.d0*z3+4.d0*z3**2-80.d0/27.d0*z5)*m2c)*amu**3
      end if

      DVA2m2 = 3.d0*sun/(4.d0*pi**2*s)

      return
      end function DVA2m2
!-----------------------------------------------------------------------
