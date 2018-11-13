!     Definition of used set of weight functions for contour integrals
      module weights_used
         implicit none

         integer, parameter :: Nmom = 1 ! Number of used weights

      contains

      function wgtr(n,x)
      use num_const
      use weights_rho

      implicit none
      integer,     intent(in) :: n
      complex(dc), intent(in) :: x
      complex(dc)             :: wgtr

      wgtr = i0

      select case (n)
      case (1)
         wgtr = wr00(x)
!        wgtr = wrp0()
!        wgtr = wrp0m1(x)
!        wgtr = wrp0mx3(x)
!        wgtr = wrp0m12(x)
      case (2)
!        wgtr = wr10(x)
!        wgtr = wrp0()
!        wgtr = wrp0m1(x)
         wgtr = wrp0m12(x)
      case (3)
         wgtr = wr11(x)
!        wgtr = wr00(x)
!        wgtr = wrp0m11(x)
!        wgtr = wrp0mx3(x)
      case (4)
!        wgtr = wr00(x)
         wgtr = wr12(x)
!        wgtr = wrp0m12(x)
!        wgtr = wrp0mx3(x)
      case (5)
!        wgtr = wr00(x)
         wgtr = wr13(x)
      case (6)
         wgtr = wrM1(x)
      end select

      end function wgtr

      function wgtD(n,x)
      use num_const
      use weights_D

      implicit none
      integer,     intent(in) :: n
      complex(dc), intent(in) :: x
      complex(dc)             :: wgtD

      wgtD = i0

      select case (n)
      case (0)
         wgtD = wDspec()
      case (1)
         wgtD = wD00(x)
!        wgtD = wDp0(x)
!        wgtD = wDp0m1(x)
!        wgtD = wDp0mx3(x)
!        wgtD = wDp0m12(x)
      case (2)
!        wgtD = wD10(x)
!        wgtD = wDp0(x)
!        wgtD = wDp0m1(x)
         wgtD = wDp0m12(x)
      case (3)
         wgtD = wD11(x)
!        wgtD = wD00(x)
!        wgtD = wDp0m11(x)
!        wgtD = wDp0mx3(x)
      case (4)
!        wgtD = wD00(x)
         wgtD = wD12(x)
!        wgtD = wDp0m12(x)
!        wgtD = wDp0mx3(x)
      case (5)
!        wgtD = wD00(x)
         wgtD = wD13(x)
      case (6)
         wgtD = wDM1(x)
      end select

      end function wgtD

      end module weights_used
