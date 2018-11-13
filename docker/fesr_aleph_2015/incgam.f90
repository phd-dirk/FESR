
!     Incomplete Gamma-function and exponential integral E
!     for complex argument z
!     Last change: Matthias Jamin, 4.6.2009

!     Adapted from Eric Kostlan & Dmitry Gokhman
!     March  1986
!     For documentation, see:
!     http://zakuski.math.utsa.edu/~gokhman/papers/igf.html

!-----------------------------------------------------------------------

      function expintE(a, z)
      use num_const

      implicit none
      real(dp),    intent(in) :: a
      complex(dc), intent(in) :: z

      complex(dc) :: expintE, incgam

      expintE = z**(a-1.d0)*incgam(1.d0-a,z)

      return
      end function expintE

!-----------------------------------------------------------------------

      function incgam(a, z)
      use num_const

      implicit none
      real(dp),    intent(in) :: a
      complex(dc), intent(in) :: z
      complex(dc)             :: incgam
 
      integer, parameter     :: ibuf = 34
      real(dp), parameter    :: zero = 0.d0, zlim = 1.d0
      complex(dc), parameter :: re = (0.36787944117144232d0, 0.d0)

      integer     :: i, ilim
      real(dp)    :: dnrm
      complex(dc) :: p, q, cdh

!     If z is near the negative real axis, then shift to z=1.
      if (dnrm(z) < zlim .or. real(z, dp) < zero .and.&
         &abs(aimag(z)) < zlim) then
         incgam = re/cdh(a, i1)
!        ilim = real(z/re, dp)
         ilim = idnint(real(z/re, dp))
         do i = 0, ibuf - ilim
            call term(a, z, i, p, q)
            incgam = incgam + p * q
         end do
      else
         incgam = exp(-z + a*log(z)) / cdh(a, z)
      end if

      return
      end function incgam

!-----------------------------------------------------------------------

      function cdh(a, z)
      use num_const

      implicit none
      real(dp),    intent(in) :: a
      complex(dc), intent(in) :: z
      complex(dc)             :: cdh

      integer     :: i, n
      real(dp)    :: rn, a1
      complex(dc) :: sun, term, cdhs

!     If Re(a-z) is too big, shift a.
!     n = real(a-z, dp)
      n = idnint(real(a-z, dp))
      if (n > 0) then
         rn = n
         a1 = a - rn
         term = i1 / z
         sun = term
         do i = 1, n - 1
            rn = n - i
            term = term * (a1 + rn) / z
            sun = term + sun
         end do 
         sun = sun + term * a1 / cdhs(a1, z)
         cdh = i1 / sun
      else
         cdh = cdhs(a, z)
      end if

      return
      end function cdh

!-----------------------------------------------------------------------

      function cdhs(a, z)
      use num_const

      implicit none
      real(dp),    intent(in) :: a
      complex(dc), intent(in) :: z
      complex(dc)             :: cdhs

      integer, parameter  :: ilim = 100000
      real(dp), parameter :: zero = 0.d0, half = 0.5d0, one = 1.d0
      real(dp), parameter :: tol1 = 1.d+10, tol2 = 1.d-10,&
                            &error = 1.d-16

      integer     :: i
      real(dp)    :: dnrm
      complex(dc) :: p0, q0, p1, q1, r0, r1, ci, factor

      q0 = one
      q1 = one
      p0 = z
      p1 = z + one - a
      do i = 1, ilim
         ci = i
         if (p0 /= zero .and. q0 /= zero .and. q1 /= zero) then
            r0 = p0 / q0
            r1 = p1 / q1
            if (dnrm(r0-r1) <= dnrm(r1)*error) then
               cdhs = r1
               return
            end if
!           Occasionally renormalise the sequences to avoid over(under)flow.
            if (dnrm(p0) > tol1 .or. dnrm(p0) < tol2 .or.&
               &dnrm(q0) > tol1 .or. dnrm(q0) < tol2) then
               factor = p0 * q0
               p0 = p0 / factor
               q0 = q0 / factor
               p1 = p1 / factor
               q1 = q1 / factor
            end if 
         end if
         p0 = z * p1 + ci * p0
         q0 = z * q1 + ci * q0
         p1 = p0 + (ci+one-a) * p1
         q1 = q0 + (ci+one-a) * q1
      end do
!     If the peripheral routines are written correctly,
!     the following four statements should never be executed.
      write(*, *) 'cdhs:  *** Warning: i >', ilim
      write(*, *) 'cdhs:  *** r0,r1= ', r0, r1
      write(*, *) ' a = ', a, '  z = ', z
      cdhs = half * (r0+r1)

      return
      end function cdhs

!-----------------------------------------------------------------------

!     Calculate p*q = -1**i(1 - x**(alpha+i))/(alpha+i)i ! carefully.
      subroutine term(a, z, i, p, q)
      use num_const

      implicit none
      real(dp),    intent(in)  :: a
      complex(dc), intent(in)  :: z
      integer,     intent(in)  :: i
      complex(dc), intent(out) :: p, q

      real(dp), parameter :: zero = 0.d0, one = 1.d0, two = 2.d0
      real(dp), parameter :: tol = 3.d-7, xlim = 39.d0

      real(dp)    :: ai, ri, dnrm
      complex(dc) :: cdlz

      if (i == 0) q = one
      ri = i
      ai = a + ri
      if (z == zero) then
         p = one / ai
         if (i /= 0) q = -q / ri
         return 
      end if
      cdlz = log(z)

!     If (1 - x**ai) = -x**ai on the computer, then
!     change the inductive scheme to avoid overflow.
      if (real(ai*cdlz, dp) > xlim .and. i /= 0) then
         p = p * (ai - one) / ai
         q = -q * z / ri
         return
      end if
      if (dnrm(dcmplx(ai,0)) > tol) then
         p = (one - z**ai) / ai
      else
         p = -cdlz * (one + cdlz*ai/two)
      end if
      if (i /= 0) q = -q / ri

      return
      end subroutine term

!-----------------------------------------------------------------------

      function dnrm(z)
      use num_const

      implicit none
      complex(dc), intent(in) :: z
      real(dp)                :: dnrm

      dnrm = abs(real(z, dp)) + abs(aimag(z))

      return
      end function dnrm

!-----------------------------------------------------------------------
