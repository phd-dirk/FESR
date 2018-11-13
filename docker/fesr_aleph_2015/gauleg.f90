
!     Routine for Gauss integration
!     Adapted from Numerical Recipes
!     Last change: Matthias Jamin, 30.10.2009

!-----------------------------------------------------------------------

      subroutine gauleg(x1,x2,x,w,n)
      use num_const

      implicit none
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x1, x2
      real(dp), intent(out) :: x(n), w(n)

      real(dp), parameter :: eps = 1.d-15

      integer  :: i, j, m
      real(dp) :: xm, xl, z, z1, p1, p2, p3, pp

      m = int(n+1)/2
      xm = 0.5d0*(x2+x1)
      xl = 0.5d0*(x2-x1)

      do i = 1,m
         z = cos(pi*(i-0.25d0)/(n+0.5d0))

   10    p1 = 1.d0
         p2 = 0.d0
         do j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         end do
         pp = n*(z*p1-p2)/(z*z-1.d0)
         z1 = z
         z = z1-p1/pp
         if (abs(z-z1) > eps) goto 10
         x(i) = xm-xl*z
         x(n+1-i) = xm+xl*z
         w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i) = w(i)
      end do

      return
      end subroutine gauleg
!-----------------------------------------------------------------------
