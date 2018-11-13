
!     Phenomenological description of the pseudoscalar
!     (ud) and (us) contribution
!     Last change: Matthias Jamin, 13.4.2010

!-----------------------------------------------------------------------

!     Pseudoscalar contribution from pion pole and exited resonances

      function deltaP(m,s0)
      use params
      use weights_used

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: s0
      real(dp)             :: deltaP

      integer :: n
      integer, parameter :: ngau = 1201
      real(dp) :: pipole, xth, spi, sun, s, rhores, breitwigner
      real(dp), dimension(ngau) :: xx, wx

      spi = mpim**2
      pipole = -4.d0*fpi**2/s0*spi/(stau+2.d0*spi)*&
               &real(wgtr(m,dcmplx(spi/s0,0)),dp)

      xth = 9.d0*spi/s0
!     Initialise Gauss integration
      call gauleg(xth,1.d0,xx,wx,ngau)

      sun = 0.d0

      do n = 1, ngau
      s = s0*xx(n)
      rhores = 2.d0/s**2*( f1p**2*m1p**4*breitwigner(s,m1p,g1p) +&
                          &f2p**2*m2p**4*breitwigner(s,m2p,g2p) )
      sun = sun - real(wgtr(m,dcmplx(xx(n),0)),dp)*2.d0*s/&
                 &(stau+2.d0*s)*rhores*wx(n)
      end do

      deltaP = 4.d0*pi**2*( pipole + sun )

      return
      end function deltaP
!-----------------------------------------------------------------------

!     Pseudoscalar contribution from Kaon pole and exited resonances

      function deltaK(m,s0)
      use params
      use weights_used

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: s0
      real(dp)             :: deltaK

      integer :: n
      integer, parameter :: ngau = 1201
      real(dp) :: kpole, xk, xth, sun, s, rhores, breitwigner
      real(dp), dimension(ngau) :: xx, wx

      xk = mkm**2/s0
      kpole = -4.d0*fk**2/s0*xk/(1.d0+2.d0*xk)*&
              &real(wgtr(m,dcmplx(xk,0)),dp)

      xth = (mkm+2.d0*mpim)**2/s0
!     Initialise Gauss integration
      call gauleg(xth,1.d0,xx,wx,ngau)

      sun = 0.d0

      do n = 1, ngau
      s = s0*xx(n)
      rhores = 2.d0/s**2*( f1k**2*m1k**4*breitwigner(s,m1k,g1k) +&
                          &f2k**2*m2k**4*breitwigner(s,m2k,g2k) )
      sun = sun - real(wgtr(m,dcmplx(xx(n),0)))*2.d0*xx(n)/&
                 &(1.d0+2.d0*xx(n))*rhores*wx(n)
      end do

      deltaK = 4.d0*pi**2*( kpole + sun )

      return
      end function deltaK
!-----------------------------------------------------------------------

!     Breit-Wigner resonance shape

      function breitwigner(s,mbw,gbw)
      use num_const

      implicit none
      real(dp), intent(in) :: s, mbw, gbw
      real(dp)             :: breitwigner

      breitwigner = mbw*gbw/pi/((s-mbw**2)**2 + mbw**2*gbw**2) 

      return
      end function breitwigner
!-----------------------------------------------------------------------
