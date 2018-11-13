
!     Defines delta_V/A^(D) as well as delta_V+A^(D)
!     Last change: Matthias Jamin, 22.1.2010

!-----------------------------------------------------------------------

!     delta_V/A^(0) in FOPT defined by taking away factor N_c/2

      function delVA0FO(m,cV0,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc)          :: delVA0FO

      complex(dc)          :: cintVA0FO

      delVA0FO = ( cintVA0FO(m,cV0,s0,mu2,astau,nmax) -&
                   cintVA0FO(m,cV0,s0,mu2,astau,0) )/1.5d0

      return
      end function delVA0FO
!-----------------------------------------------------------------------

!     delta_V/A^(0) in CIPT defined by taking away factor N_c/2

      function delVA0CI(m,cV0,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc)          :: delVA0CI

      complex(dc)          :: cintVA0CI

      delVA0CI = ( cintVA0CI(m,cV0,s0,mu2,astau,nmax) -&
                   cintVA0CI(m,cV0,s0,mu2,astau,0) )/1.5d0

      return
      end function delVA0CI
!-----------------------------------------------------------------------

!     delta_V/A^(0) in BS defined by taking away factor N_c/2

      function delVA0BS(m,c51,s0,astau)
      use num_const

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: c51, s0, astau
      complex(dc)          :: delVA0BS

      complex(dc)          :: cintVA0BS

      delVA0BS = cintVA0BS(m,c51,s0,astau)/1.5d0

      return
      end function delVA0BS
!-----------------------------------------------------------------------

!     delta_V+A^(0) in FOPT defined by taking away factor N_c

      function delVpA0FO(m,cV0,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc)          :: delVpA0FO

      complex(dc)          :: cintVpA0FO

      delVpA0FO = ( cintVpA0FO(m,cV0,s0,mu2,astau,nmax) -&
                    cintVpA0FO(m,cV0,s0,mu2,astau,0) )/3.d0

      return
      end function delVpA0FO
!-----------------------------------------------------------------------

!     delta_V+A^(0) in CIPT defined by taking away factor N_c

      function delVpA0CI(m,cV0,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, nmax
      real(dp), intent(in) :: s0, mu2, astau
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      complex(dc)          :: delVpA0CI

      complex(dc)          :: cintVpA0CI

      delVpA0CI = ( cintVpA0CI(m,cV0,s0,mu2,astau,nmax) -&
                    cintVpA0CI(m,cV0,s0,mu2,astau,0) )/3.d0

      return
      end function delVpA0CI
!-----------------------------------------------------------------------

!     delta_V+A^(0) in BS defined by taking away factor N_c

      function delVpA0BS(m,c51,s0,astau)
      use num_const

      implicit none
      integer,  intent(in) :: m
      real(dp), intent(in) :: c51, s0, astau
      complex(dc)          :: delVpA0BS

      complex(dc)          :: cintVpA0BS

      delVpA0BS = cintVpA0BS(m,c51,s0,astau)/3.d0

      return
      end function delVpA0BS
!-----------------------------------------------------------------------

!     delta_V/A^(2) in FOPT defined by taking away factor N_c/2

      function delVA2FO(m,r,i,j,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, r, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: delVA2FO

      complex(dc)          :: cintVA2FO

      delVA2FO = cintVA2FO(m,r,i,j,s0,mu2,astau,nmax)/1.5d0

      return
      end function delVA2FO
!-----------------------------------------------------------------------

!     delta_V+A^(2) in FOPT defined by taking away factor N_c

      function delVpA2FO(m,i,j,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: delVpA2FO

      complex(dc)          :: cintVpA2FO

      delVpA2FO = cintVpA2FO(m,i,j,s0,mu2,astau,nmax)/3.d0

      return
      end function delVpA2FO
!-----------------------------------------------------------------------

!     delta_V/A^(2) in CIPT defined by taking away factor N_c/2

      function delVA2CI(m,r,i,j,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, r, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: delVA2CI

      complex(dc)          :: cintVA2CI

      delVA2CI = cintVA2CI(m,r,i,j,s0,mu2,astau,nmax)/1.5d0

      return
      end function delVA2CI
!-----------------------------------------------------------------------

!     delta_V+A^(2) in CIPT defined by taking away factor N_c

      function delVpA2CI(m,i,j,s0,mu2,astau,nmax)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j, nmax
      real(dp), intent(in) :: s0, mu2, astau
      complex(dc)          :: delVpA2CI

      complex(dc)          :: cintVpA2CI

      delVpA2CI = cintVpA2CI(m,i,j,s0,mu2,astau,nmax)/3.d0

      return
      end function delVpA2CI
!-----------------------------------------------------------------------

!     delta_V/A^(4) in FOPT defined by taking away factor N_c/2

      function delVA4FO(m,r,i,j,s0,mu2,astau,aGGinv)
      use num_const

      implicit none
      integer,  intent(in) :: m, r, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: delVA4FO

      complex(dc)          :: cintVA4FO

      delVA4FO = cintVA4FO(m,r,i,j,s0,mu2,astau,aGGinv)/1.5d0

      return
      end function delVA4FO
!-----------------------------------------------------------------------

!     delta_V+A^(4) in FOPT defined by taking away factor N_c

      function delVpA4FO(m,i,j,s0,mu2,astau,aGGinv)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: delVpA4FO

      complex(dc)          :: cintVpA4FO

      delVpA4FO = cintVpA4FO(m,i,j,s0,mu2,astau,aGGinv)/3.d0

      return
      end function delVpA4FO
!-----------------------------------------------------------------------

!     delta_V/A^(4) in CIPT defined by taking away factor N_c/2

      function delVA4CI(m,r,i,j,s0,mu2,astau,aGGinv)
      use num_const

      implicit none
      integer,  intent(in) :: m, r, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: delVA4CI

      complex(dc)          :: cintVA4CI

      delVA4CI = cintVA4CI(m,r,i,j,s0,mu2,astau,aGGinv)/1.5d0

      return
      end function delVA4CI
!-----------------------------------------------------------------------

!     delta_V+A^(4) in CIPT defined by taking away factor N_c

      function delVpA4CI(m,i,j,s0,mu2,astau,aGGinv)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j
      real(dp), intent(in) :: s0, mu2, astau, aGGinv
      complex(dc)          :: delVpA4CI

      complex(dc)          :: cintVpA4CI

      delVpA4CI = cintVpA4CI(m,i,j,s0,mu2,astau,aGGinv)/3.d0

      return
      end function delVpA4CI
!-----------------------------------------------------------------------

!     delta_V/A^(6,8) defined by taking away factor N_c/2

      function delVA68(m,r,i,j,s0,astau,rhoD6,c8D8)
      use num_const

      implicit none
      integer,  intent(in) :: m, r, i, j
      real(dp), intent(in) :: s0, astau, rhoD6, c8D8
      complex(dc)          :: delVA68

      complex(dc)          :: cintVA68

      delVA68 = cintVA68(m,r,i,j,s0,astau,rhoD6,c8D8)/1.5d0

      return
      end function delVA68
!-----------------------------------------------------------------------

!     delta_V+A^(6,8) defined by taking away factor N_c

      function delVpA68(m,i,j,s0,astau,rhoD6,c8D8)
      use num_const

      implicit none
      integer,  intent(in) :: m, i, j
      real(dp), intent(in) :: s0, astau, rhoD6, c8D8
      complex(dc)          :: delVpA68

      complex(dc)          :: cintVpA68

      delVpA68= cintVpA68(m,i,j,s0,astau,rhoD6,c8D8)/3.d0

      return
      end function delVpA68
!-----------------------------------------------------------------------
