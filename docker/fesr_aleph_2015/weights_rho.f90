!     Definition of weight functions for contour integrals of rho(s)
      module weights_rho
         implicit none

      contains

!     Polynomial weight (0,0)
         function wr00(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr00

         wr00 = (i1-x)**2*(i1+2.d0*x)
         end function wr00

!     Polynomial weight (0,1)
         function wr01(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr01

         wr01 = (i1-x)**2*(i1+2.d0*x)*x
         end function wr01

!     Polynomial weight (0,2)
         function wr02(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr02

         wr02 = (i1-x)**2*(i1+2.d0*x)*x**2
         end function wr02

!     Polynomial weight (0,3)
         function wr03(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr03

         wr03 = (i1-x)**2*(i1+2.d0*x)*x**3
         end function wr03

!     Polynomial weight (1,0)
         function wr10(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr10

         wr10 = (i1-x)**3*(i1+2.d0*x)
         end function wr10

!     Polynomial weight (1,1)
         function wr11(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr11

         wr11 = (i1-x)**3*(i1+2.d0*x)*x
         end function wr11

!     Polynomial weight (1,2)
         function wr12(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr12

         wr12 = (i1-x)**3*(i1+2.d0*x)*x**2
         end function wr12

!     Polynomial weight (1,3)
         function wr13(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr13

         wr13 = (i1-x)**3*(i1+2.d0*x)*x**3
         end function wr13

!     Polynomial weight (2,0)
         function wr20(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr20

         wr20 = (i1-x)**4*(i1+2.d0*x)
         end function wr20

!     Polynomial weight (2,1)
         function wr21(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr21

         wr21 = (i1-x)**4*(i1+2.d0*x)*x
         end function wr21

!     Polynomial weight (2,2)
         function wr22(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr22

         wr22 = (i1-x)**4*(i1+2.d0*x)*x**2
         end function wr22

!     Polynomial weight (2,3)
         function wr23(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wr23

         wr23 = (i1-x)**4*(i1+2.d0*x)*x**3
         end function wr23

!     Power weight (0)
         function wrp0()
         use num_const

         implicit none
         complex(dc)             :: wrp0

         wrp0 = i1
         end function wrp0

!     Power weight (1)
         function wrp1(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp1

         wrp1 = x
         end function wrp1

!     Power weight (0-1)
         function wrp0m1(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m1

         wrp0m1 = i1-x
         end function wrp0m1

!     Power weight (0-2*1)
         function wrp0m21(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m21

         wrp0m21 = i1-2.d0*x
         end function wrp0m21

!     Power weight (2)
         function wrp2(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp2

         wrp2 = x*x
         end function wrp2

!     Weight function 1-x^2
         function wrp0m11(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m11

         wrp0m11 = i1-x**2
         end function wrp0m11

!     Weight function (1-x)^2
         function wrp0m12(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m12

         wrp0m12 = (i1-x)**2
         end function wrp0m12

!     Weight function 1-x^3
         function wrp0mx3(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0mx3

         wrp0mx3 = i1-x**3
         end function wrp0mx3

!     Weight function (1-x)^3
         function wrp0m13(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m13

         wrp0m13 = (i1-x)**3
         end function wrp0m13

!     Weight function (1-x)^3 (1+3x)
         function wrp0m13nx(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m13nx

         wrp0m13nx = (i1-x)**3*(i1+3.d0*x)
         end function wrp0m13nx

!     Weight function (1-x)^4 (1+4x)
         function wrp0m14nx(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m14nx

         wrp0m14nx = (i1-x)**4*(i1+4.d0*x)
         end function wrp0m14nx

!     Weight function (1-x)^2 exp(-x)
         function wrp0m12exp(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrp0m12exp

         wrp0m12exp = (i1-x)**2/exp(x)
         end function wrp0m12exp

!     Maltman weight w_0(x)
         function wrM0(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrM0

         wrM0 = x*(i1-x)**2
         end function wrM0

!     Maltman weight w_N(x)
         function wrMN(n,x)
         use num_const

         implicit none
         integer, intent(in) :: n
         complex(dc), intent(in) :: x
         complex(dc)             :: wrMN

         wrMN = i1 - n/(n-i1)*x + x**n/(n-i1)
         end function wrMN

!     Martin's first weight
         function wrM1(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrM1

         wrM1 = 2.d0/3.d0*((i1-x)**2 + 2.d0/pi*sin(pi*x))

         end function wrM1

!     Martin's second weight
         function wrM2(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrM2

         wrM2 = 2.d0/3.d0*(i1+x**2-pi**2*x**3/3.d0+(-6.d0+&
         &2.d0*pi**2/3.d0)*x**5+(4.d0-pi**2/3.d0)*x**7)

         end function wrM2

!     Martin's third weight
         function wrM3(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wrM3

         wrM3 = 2.d0/3.d0*(i1+2.d0*x**2-4.d0*x**3-3.d0*x**4+4.d0*x**5)

         end function wrM3

      end module weights_rho

