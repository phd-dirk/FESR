!     Useful matrix manipulation routines
      module matrix
         implicit none

      contains

!-----------------------------------------------------------------------

!     Subroutine to find the inverse of a square matrix
!     Author : Louisda16th a.k.a Ashwith J. Rego
!     Reference : Algorithm has been well explained in:
!     http://math.uww.edu/~mcfarlat/inverse.htm           
!     http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

      subroutine matinv(matrix, inverse, n, errorflag)
      use num_const

      implicit none
      integer, intent(in)  :: n
      integer, intent(out) :: errorflag  ! Return error status. -1 for error, 0 for normal
      real(dp), intent(in),  dimension(n,n) :: matrix  ! Input matrix
      real(dp), intent(out), dimension(n,n) :: inverse ! Inverted matrix

      logical :: flag = .true.

      integer :: i, j, k
      real(dp) :: m
      real(dp), dimension(n,2*n) :: augmatrix ! augmented matrix
      
!     Augment input matrix with an identity matrix
      do i = 1, n
         do j = 1, 2*n
            if (j <= n ) then
               augmatrix(i,j) = matrix(i,j)
            else if ((i+n) == j) then
               augmatrix(i,j) = 1.d0
            else
               augmatrix(i,j) = 0.d0
            end if
         end do
      end do

!     Reduce augmented matrix to upper triangular form
      do k =1, n-1
         if (augmatrix(k,k) == 0.d0) then
            flag = .false.
            do i = k+1, n
               if (augmatrix(i,k) /= 0.d0) then
                  do j = 1, 2*n
                     augmatrix(k,j) = augmatrix(k,j) + augmatrix(i,j)
                  end do
                  flag = .true.
                  exit
               end if
               if (flag .eqv. .false.) then
                  print*, "Matrix is non - invertible"
                  inverse = 0.d0
                  errorflag = -1
                  return
               end if
            end do
         end if
         do j = k+1, n
            m = augmatrix(j,k)/augmatrix(k,k)
            do i = k, 2*n
               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            end do
         end do
      end do

!     Test for invertibility
      do i = 1, n
         if (augmatrix(i,i) == 0.d0) then
            print*, "Matrix is non - invertible"
            inverse = 0.d0
            errorflag = -1
            return
         end if
      end do

!     Make diagonal elements as 1
      do i = 1, n
         m = augmatrix(i,i)
         do j = i, (2*n)
            augmatrix(i,j) = (augmatrix(i,j) / m)
         end do
      end do

!     Reduced right side half of augmented matrix to identity matrix
      do k = n-1, 1, -1
        do i =1, k
           m = augmatrix(i,k+1)
           do j = k, (2*n)
              augmatrix(i,j) = augmatrix(i,j) - augmatrix(k+1,j) * m
           end do
         end do
      end do

!     Store answer
      do i =1, n
         do j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
         end do
      end do
      errorflag = 0

      end subroutine matinv

!-----------------------------------------------------------------------

!     Given a positive definite matrix a(1:n,1:n), with physical dimension  
!     np, this routine constructs its Cholesky decomposition, A = LL^T.  On  
!     input, only the upper triangle of A need be given; it is not modified.
!     The Cholesky factor L is returned in the lower triangle of A, except  
!     for its diagonal elements, which are returned in p(1:n).              
!     (Numerical Recipes, Second Edition, Sec. 2.9)                         

      subroutine choldc(a,n,np,p) 
      use num_const

      implicit none 
      integer,  intent(in) :: n, np 
      real(dp)             :: a(np,np), p(n)

      integer  :: i, j, k 
      real(dp) :: sun 
 
      do i = 1, n 
         do j = i, n 
            sun = a(i,j) 
            do k = i-1, 1, -1 
               sun = sun - a(i,k) * a(j,k) 
            end do 
            if (i.eq.j) then 
               if (sun <= 0) write(*,*) 'choldc failed' 
               p(i) = sqrt(sun) 
            else 
               a(j,i) = sun/p(i) 
            end if 
         end do 
      end do 
      return 

      end subroutine choldc

!-----------------------------------------------------------------------

!     Solves the set of n linear equations Lx = b, with L a lower triangular
!     matrix.  b(1:n) is the input as the right-hand side vector b, and     
!     returns the solution vector x. L, n and np are not modified by this   
!     routine, and can be left in place for successive calls with different 
!     right-hand sides b.                                                   
!     (Adapted from Numerical Recipes, Second Edition, Sec. 2.3)            

      subroutine lbksb(L,n,np,b) 
      use num_const

      integer,  intent(in) :: n, np 
      real(dp), intent(in) :: L(np,np)
      real(dp)             :: b(n) 

      integer i, j 
      real(dp) sun 

      b(1) = b(1)/L(1,1) 
      do i = 2, n 
         sun = b(i) 
         do j = 1, i-1 
            sun = sun - L(i,j)*b(j) 
         end do 
         b(i) = sun/L(i,i) 
      end do 
      return 

      end subroutine lbksb

!-----------------------------------------------------------------------

!     Test matrix inverse

      subroutine matinvtest(matrix,matinv,n)
      use num_const

      implicit none 
      integer  :: n, i, j, k
      real(dp) :: matrix(n,n), matinv(n,n), onemat(n,n)

      do i = 1, n
         do j = 1, n
            onemat(i,j) = 0.d0
         end do
      end do

      do i = 1, n
         do j = 1, n
            do k = 1, n
               onemat(i,j) = onemat(i,j) + matinv(i,k)*matrix(k,j)
            end do
         end do
      end do

      write(*,*) ""
      do i = 1, n
         write(*,*) (/ (onemat(i,j), j=1,n) /)
      end do

      return
      end subroutine matinvtest
!-----------------------------------------------------------------------

      subroutine milc_proc(n,nevs,corin,corout)
      use num_const

      implicit none
      integer,  intent(in) :: n, nevs
      real(dp), intent(in) :: corin(n,n)

      integer  :: i, j, k, ierr
      real(dp) :: corout(n,n), eival(n), sdiag(n), eivec(n,n), sun

      call tred2(n,n,corin,eival,sdiag,eivec)
      call tql2(n,n,eival,sdiag,eivec,ierr)

      if (ierr /= 0) write(8,*) 'Error in eigen decomposition !!!'

!     print eigenvalues
      write(8,*) 'Eigenvalues of correlation matrix:'
      do i = 1, n
         write(8,*) i, eival(i)
      end do
      write(8,*) ''

!     Set eigenvalues beyond nevs to zero
      do i = n-nevs, 1, -1
         eival(i) = 0.d0
      end do

      do i = 1, n
         do j = 1, n
            sun = 0.d0
            do k = 1, n
               sun = sun + eivec(i,k)*eival(k)*eivec(j,k) 
            end do
            corout(i,j) = sun
         end do
      end do

      do i = 1, n
         corout(i,i) = 1.d0
      end do

      return
      end subroutine milc_proc
!-----------------------------------------------------------------------

      subroutine tql2(nm,n,d,e,z,ier)
      use num_const

!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO tred2.
!
!     CALLING MODE:
!               CALL tql2(nm,n,d,e,z,ier)
!     INPUTS:
!     nm  (I4)  1ST DIMENSION OF MATRICES a AND z IN CALLING PROGRAM
!     n   (I4)  SIZE OF z
!     d  (R*8)  MAIN DIAGONAL (n) OF THE TRIDIAGONAL MATRIX
!     e  (R*8)  SUB-DIAGONAL  (n) OF THE TRIDIAGONAL MATRIX
!     z  (R*8)  TABLE (nm,n) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY d AND e, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO tred2.
!     OUTPUTS:
!     d  (R*8)  EIGENVALUES
!     z  (R*8)  EIGENVECTORS
!     ier (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.

      implicit none 
      integer  :: i, j, k, l, m, n, nm, jm, ier
      real(dp) :: d(n), e(n), z(nm,n), b, c, f, g, h, p, r, s, eps, eps1

      data eps /0.d0/, jm /30/
      ier = 0
      if (n == 1) go to 38
!
!     MACHINE EPSILON
!
      if (eps /= 0.d0) go to 12
      eps = 1.d0
   10 eps = eps/2.d0
      eps1 = 1.d0 + eps
      if (eps1 > 1.d0) go to 10
!
   12 do 14 i = 2, n
   14 e(i-1) = e(i)
      e(n) = 0.d0
      f = 0.d0
      b = 0.d0
!
      do 28 l = 1, n
      j = 0
      h = eps*(abs(d(l)) + abs(e(l)))
      if (b < h) b = h
!
!     SEEK SMALLEST ELEMENT OF SUBDIAGONAL
!
      do 16 m = l, n
      if (abs(e(m)) <= B) go to 18
   16 continue
   18 if (m == l) go to 26

!     START ITERATION

   20 if (j == jm) go to 36
      j = j+1

!     SHIFT

      g = d(l)
      p = (d(l+1) - g)/(2.d0*e(l))
      r = sqrt(p*p + 1.d0)
      d(l) = e(l)/(p + sign(r,p))
      h = g - d(l)
      do 22 i = l+1, n
   22 d(i) = d(i) - h
      f = f + h

!     QL TRANSFORMATION

      p = d(m)
      c = 1.d0
      s = 0.d0
      do 24 i = m-1, l, -1
      g = c*e(I)
      h = c*p
      if (abs(p) >= abs(e(i))) then
         c = e(i)/p
         r = sqrt(c*c + 1.d0)
         e(i+1) = s*p*r
         s = c/r
         c = 1.d0/r
      else
         c = p/e(i)
         r = sqrt(c*c + 1.d0)
         e(i+1) = s*e(i)*r
         s = 1.d0/r
         c = c*s
      end if
      p = c*d(i) - s*g
      d(i+1) = h + s*(c*g + s*d(i))

!     ELEMENTS OF EIGENVECTORS

      do 24 k = 1, n
      h = z(k,i+1)
      z(k,i+1) = s*z(k,i) + c*h
      z(k,i) = z(k,i)*c - s*h
   24 continue
      e(L) = s*p
      d(L) = c*p
      if (abs(e(l)) > b) go to 20

!     CONVERGENCE

   26 d(l) = d(l) + f
   28 continue

!     SORT EIGENVALUES AND EIGENVECTORS
!     IN ASCENDING ORDER

      do 34 l = 2,n
      i = l-1
      k = i
      p = d(i)
      do 30 j = l, n
      if (d(j) >= p) go to 30
      k = j
      p = d(j)
   30 continue
      if (k == i) go to 34
      d(k) = d(i)
      d(i) = p
      do 32 j = 1, n
      p = z(j,i)
      z(j,i) = z(j,k)
   32 z(j,k) = p
   34 continue
      go to 38

!     NO CONVERGENCE

   36 ier = l
   38 return
      end subroutine tql2
!-----------------------------------------------------------------------

      subroutine tred2(nm,n,a,d,e,z)
      use num_const

!     TRIDIAGONALIZATION OF A SYMMETRIC MATRIX BY ORTHOGONAL TRANSFORMATIONS
!     (ALGORITHM OF HOUSEHOLDER)
!     CALLING MODE:
!               CALL tred2(nm,n,a,d,e,z)
!     INPUTS:
!     nm  (I4)  1ST DIMENSION OF MATRICES a AND a IN CALLING PROGRAM
!     n   (I4)  SIZE OF a
!     a  (R*8)  TABLE(nm,n) STORING THE COEFFICIENTS OF SYMMETRIC a MATRIX
!               (LOWER HALF), a IS NOT DESTROYED DURING THE PROCESS
!               IF z MATRIX HAS NOT THE SAME ADDRESS.
!     OUTPUTS:
!     d  (R*8)  MAIN DIAGONAL (n) OF REDUCED TRIDIAGONAL MATRIX
!     e  (R*8)  SUB-DIAGONAL  (n) OF REDUCED TRIDIAGONAL MATRIX
!     z  (R*8)  TABLE (nm,n) STORING THE ELEMENTS OF THE ORTHOGONAL 
!               TRANSFORMATION MATRIX.
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.

      implicit none 
      integer  :: i, j, k, l, n, nm
      real(dp) :: a(nm,n), d(n), e(n), z(nm,n), f, g, h, hh, scale

!     LOWER HALF OF a PUT INTO z

      do 10 i = 1, n
      do 10 j = 1, i
   10 z(i,j) = a(i,j)
      if (n == 1) go to 32

!     n-2 STAGE OF TRANSFORMATION

      do 30 i = n, 2, -1
      l = i-1
      h = 0.d0

!     CONDITIONNING BY NORM OF a

      scale = 0.d0
      if (l < 2) go to 14
      do 12 k = 1, l
   12 scale = scale + abs(z(i,k))
      if (scale /= 0.d0) go to 16

   14 e(i) = z(i,l)
      go to 28

   16 do 18 k = 1, l
      z(i,k) = z(i,k)/scale
      h = h + z(i,k)*z(i,k)
   18 continue

      f = z(i,l)
      g = -sign(sqrt(h),f)
      e(i) = scale*g
      h = h - f*g
      z(i,l) = f - g
      f = 0.d0
      do 24 j = 1, l
      z(j,i) = z(i,j)/h
      g = 0.d0

!     ELEMENT OF a*u
      do 20 k = 1, j
   20 g = g + z(j,k)*z(i,k)
      if (l >= j+1) then
      do 22 k = j+1, l
   22 g = g + z(k,j)*z(i,k)

!     ELEMENT OF p = a*u/h

      end if
      e(j) = g/h
      f = f + e(j)*z(i,j)
   24 continue

!     ELEMENT OF k

      hh = f/(h + h)

!     REDUCED FORM OF a

      do 26 j = 1, l
      f = z(i,j)
      g = e(j) - hh*f
      e(j) = g
      do 26 k = 1, j
      z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
   26 continue
!
   28 d(i) = h
   30 continue

!     END OF TRANSFORMATION

   32 d(1) = 0.d0
      e(1) = 0.d0

!     ACCUMULATE TRANSFORMATION MATRICES IN z

      DO 40 i = 1, n
      l = i-1
      if (d(i) /= 0.d0) then
      do 36 j = 1, l
      g = 0.d0
      do 34 k = 1, l
   34 g = g + z(i,k)*z(k,j)
      do 36 k = 1, l
      z(k,j) = z(k,j) - g*z(k,i)
   36 continue
      end if
      d(i) = z(i,i)
      z(i,i) = 1.d0
      if (l < 1) go to 40
      do 38 j = 1, l
      z(i,j) = 0.d0
      z(j,i) = 0.d0
   38 continue
   40 continue

      return
      end subroutine tred2
!-----------------------------------------------------------------------
      end module matrix
