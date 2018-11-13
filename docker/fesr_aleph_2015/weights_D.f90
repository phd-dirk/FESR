!     Definition of weight functions for contour integrals of D(s)
      module weights_D
         implicit none

      contains

!     Weight to yield spectral function
         function wDspec()
         use num_const

         implicit none
         complex(dc)             :: wDspec

         wDspec = i1
         end function wDspec

!     Polynomial weight (0,0)
         function wD00(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD00

         wD00 = (i1-x)**3*(i1+x)
         end function wD00

!     Polynomial weight (0,1)
         function wD01(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD01

         wD01 = (i1-x)**3*(3.d0+9.d0*x+8.d0*x**2)/1.d1
         end function wD01

!     Polynomial weight (0,2)
         function wD02(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD02

         wD02 = 2.d0*(i1-x)**3*(1.d0+3.d0*x+6.d0*x**2+5.d0*x**3)/15.d0
         end function wD02

!     Polynomial weight (0,3)
         function wD03(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD03

         wD03 = (i1-x)**3*(1.d0+3.d0*x+6.d0*x**2+1.d1*x**3+&
                          &8.d0*x**4)/14.d0
         end function wD03

!     Polynomial weight (1,0)
         function wD10(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD10

         wD10 = (i1-x)**4*(7.d0+8.d0*x)/1.d1
         end function wD10

!     Polynomial weight (1,1)
         function wD11(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD11

         wD11 = (i1-x)**4*(1.d0+2.d0*x)**2/6.d0
         end function wD11

!     Polynomial weight (1,2)
         function wD12(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD12

         wD12 = (i1-x)**4*(13.d0+52.d0*x+13.d1*x**2+12.d1*x**3)/21.d1
         end function wD12

!     Polynomial weight (1,3)
         function wD13(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD13

         wD13 = (i1-x)**4*(2.d0+8.d0*x+2.d1*x**2+4.d1*x**3+&
                          &35.d0*x**4)/7.d1
         end function wD13

!     Polynomial weight (2,0)
         function wD20(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD20

         wD20 = 2.d0*(i1-x)**5*(4.d0+5.d0*x)/15.d0
         end function wD20

!     Polynomial weight (2,1)
         function wD21(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD21

         wD21 = (i1-x)**5*(11.d0+55.d0*x+6.d1*x**2)/105.d0
         end function wD21

!     Polynomial weight (2,2)
         function wD22(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD22

         wD22 = (i1-x)**5*(1.d0+5.d0*x+15.d0*x**2+15.d0*x**3)/3.d1
         end function wD22

!     Polynomial weight (2,3)
         function wD23(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wD23

         wD23 = (i1-x)**5*(17.d0+85.d0*x+255.d0*x**2+595.d0*x**3+&
                         &56.d1*x**4)/126.d1
         end function wD23

!     Power weight (0)
         function wDp0(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0

         wDp0 = 2.d0*(i1-x)
         end function wDp0

!     Power weight (1)
         function wDp1(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp1

         wDp1 = i1-x*x
         end function wDp1

!     Power weight (0-1)
         function wDp0m1(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m1

         wDp0m1 = (i1-x)**2
         end function wDp0m1

!     Power weight (0-2*1)
         function wDp0m21(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m21

         wDp0m21 = 2.d0*x*(x-i1)
         end function wDp0m21

!     Power weight (2)
         function wDp2(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp2

         wDp2 = 2.d0/3.d0*(i1-x**3)
         end function wDp2

!     Weight function 1-x^2
         function wDp0m11(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m11

         wDp0m11 = 2.d0/3.d0*(i1-x)**2*(2.d0+x)
         end function wDp0m11

!     Weight function (1-x)^2
         function wDp0m12(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m12

         wDp0m12 = 2.d0/3.d0*(i1-x)**3
         end function wDp0m12

!     Weight function 1-x^3
         function wDp0mx3(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0mx3

         wDp0mx3 = (i1-x)**2*(3.d0+2.d0*x+x*x)/2.d0
         end function wDp0mx3

!     Weight function (1-x)^3
         function wDp0m13(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m13

         wDp0m13 = (i1-x)**4/2.d0
         end function wDp0m13

!     Weight function (1-x)^3 (1+3x)
         function wDp0m13nx(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m13nx

         wDp0m13nx = 2.d0/5.d0*(i1-x)**4*(2.d0+3.d0*x)
         end function wDp0m13nx

!     Weight function (1-x)^4 (1+4x)
         function wDp0m14nx(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m14nx

         wDp0m14nx = 2.d0/3.d0*(i1-x)**5*(1.d0+2.d0*x)
         end function wDp0m14nx

!     Weight function (1-x)^2 exp(-x)
         function wDp0m12exp(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDp0m12exp

         wDp0m12exp = 2.d0*(-2.d0/exp(0.d0)+(i1+x*x)/exp(x))
         end function wDp0m12exp

!     Maltman weight w_0(x)
         function wDM0(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDM0

         wDM0 = (i1-x)**3*(i1+3.d0*x)/6.d0
         end function wDM0

!     Maltman weight w_N(x)
         function wDMN(n,x)
         use num_const

         implicit none
         integer, intent(in) :: n
         complex(dc), intent(in) :: x
         complex(dc)             :: wDMN

         wDMN = n/(n+i1)-2.d0*x+n/(n-i1)*x*x-2.d0/(n*n-i1)*x**(n+1) 
         end function wDMN

!     Martin's first weight
         function wDM1(x)
         use num_const

         implicit none
         complex(dc), intent(in) :: x
         complex(dc)             :: wDM1

         wDM1 = 4.d0/9.d0*(i1-x)**3 + 8.d0/(3.d0*pi**2)*(1.d0+cos(pi*x))

         end function wDM1

      end module weights_D
