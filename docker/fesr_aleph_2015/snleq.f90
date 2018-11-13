!
!     Based on   J.J. More  and  M.Y. Cosnard
!
!       ALGORITHM 554 BRENTM, A Fortran Subroutine for the
!       Numerical Solution of Systems of Nonlinear Equations [C5]
!
!     ACM Trans. Math. Software 6 (1980) 240-251.
!
!     Adapted from CERNlib
!     Last change: Matthias Jamin, 2.3.2009

      subroutine my_snleq(n,x,f,ftol,xtol,maxf,iprt,info,sub,w)
      use num_const
      implicit none

      logical :: lcv
      integer :: n, maxf, iprt, info, i, j, k, m, mopt, mpt(288),&
                &iflag, numf, nfcall, nier6, nier7, nier8, nsing
      real(dp) :: x(n), f(n), ftol, xtol, w(n,*)
      real(dp) :: delta, difit, difit1, eps, eta, fky, fkz,&
                 &fnorm, fnorm1, h, p05, scale, sknorm, temp,&
                 &xnorm, z1

      parameter (z1 = 1.d0, scale = 1.d1, p05 = 5*z1/1.d2)
!**** eps = sqrt(smallest FP number)
!     eps on gfortran      
      parameter (eps = 0.105367121277235087d-07)
      data (mpt(i), i=1,288)&
     &/ 1*1, 1*2, 3*3, 3*4, 4*5, 4*6, 4*7, 4*8, 5*9,5*10,5*11,5*12,&
      &5*13,5*14,6*15,6*16,5*17,6*18,6*19,6*20,7*21,6*22,6*23,7*24,&
      &6*25,7*26,6*27,7*28,7*29,7*30,7*31,7*32,7*33,7*34,7*35,7*36,&
      &8*37,7*38,7*39,8*40,7*41,8*42,7*43,8*44,8*45,7*46,8*47,8*48/

      info = 0
      if(n <= 0 .or. ftol <= 0.d0 .or. xtol <= 0.d0) return
!
!     Find optimal 'mopt' for iterative refinement
!
      mopt = 1
      if (n <= 288) then
         mopt = mpt(n)
      else
         h = 0
         do 1 i = 49,n
            temp = log(i+z1)/dble(n+2*i+1)
            if(temp < h) then
               mopt = i-1
               goto 2
            end if
    1    h = temp
      end if

    2 iflag = 0
      numf = 0
      nfcall = 0

      nier6 = -1
      nier7 = -1
      nier8 = 0
      fnorm = 0.d0
      difit = 0.d0
      xnorm = 0.d0
      do i = 1,n
         xnorm = max(xnorm,abs(x(i)))
      end do
      delta = scale*xnorm
      if (xnorm == 0.d0) delta = scale

   20 if (iprt /= 0) write(6,'(1X,I5,D25.14)') (i,x(i), i=1,n)

      nsing = n
      fnorm1 = fnorm
      difit1 = difit
      fnorm = 0.d0
!
!     Compute step h for the divided difference which approximates
!     the k-th row of the Jacobian matrix
!
      h = eps*xnorm
      if (h == 0.d0) h = eps
      do j = 1,n
         do i = 1,n
            w(i,j+3) = 0.d0
         end do
         w(j,j+3) = h
         w(j,2) = x(j)
      end do
!
!     Enter a subiteration
!
      do 150 k = 1,n
      iflag = k
      call sub(n,w(1,2),f,iflag)
      fky = f(k)
      nfcall = nfcall+1
      numf = nfcall/n
      if (iflag < 0) goto 230
      fnorm = max(fnorm,abs(fky))
!
!     Compute the k-th row of the Jacobian matrix
!
      do j = k,n
         do i = 1,n
            w(i,3) = w(i,2)+w(i,j+3)
         end do
         call sub(n,w(1,3),f,iflag)
         fkz = f(k)
         nfcall = nfcall+1
         numf = nfcall/n
         if (iflag < 0) goto 230
         w(j,1) = fkz-fky
      end do
      f(k) = fky
!
!     Compute the Householder transformation to reduce the k-th row
!     of the Jacobian matrix to a multiple of the k-th unit vector
!
      eta = 0.d0
      do i = k,n
         eta = max(eta,abs(w(i,1)))
      end do
      if(eta == 0) goto 150
      nsing = nsing-1
      sknorm = 0.d0
      do i = k,n
         w(i,1) = w(i,1)/eta
         sknorm = sknorm+w(i,1)**2
      end do
      sknorm = sqrt(sknorm)
      if(w(k,1) < 0.d0) sknorm = -sknorm
      w(k,1) = w(k,1)+sknorm
!
!     Apply the transformation
!
      do i = 1,n
         w(i,3) = 0.d0
      end do
      do j = k,n
         do i = 1,n
            w(i,3) = w(i,3)+w(j,1)*w(i,j+3)
         end do
      end do
      do j = k,n
         temp = w(j,1)/(sknorm*w(k,1))
         do i = 1,n
            w(i,j+3) = w(i,j+3)-temp*w(i,3)
         end do
      end do
!
!     Compute the subiterate
!
      w(k,1) = sknorm*eta
      temp = fky/w(k,1)
      if (h*abs(temp) > delta) temp = sign(delta/h,temp)
      do i = 1,n
         w(i,2) = w(i,2)+temp*w(i,k+3)
      end do
  150 continue
!
!     Compute the norms of the iterate and correction vector
!
      xnorm = 0.d0
      difit = 0.d0
      do i = 1,n
         xnorm = max(xnorm,abs(w(i,2)))
         difit = max(difit,abs(x(i)-w(i,2)))
         x(i) = w(i,2)
      end do
!
!     Update the bound on the correction vector
!
      delta = max(delta,scale*xnorm)
!
!     Determine the progress of the iteration
!
      lcv = (fnorm < fnorm1).and.(difit < difit1).and.(nsing == 0)
      nier6 = nier6+1
      nier7 = nier7+1
      nier8 = nier8+1
      if(lcv) nier6 = 0
      if (fnorm < fnorm1 .or. difit < difit1) nier7 = 0
      if (difit > eps*xnorm) nier8 = 0
!
!     Tests for convergence
!
      if (fnorm <= ftol) info = 1
      if (difit <= xtol*xnorm .and. lcv) info = 2
      if (fnorm <= ftol .and. info == 2) info = 3
      if (info /= 0) goto 230
!
!     Tests for termination
!
      if (numf >= maxf) info = 4
      if (nsing == n) info = 5
      if (nier6 == 5) info = 6
      if (nier7 == 3) info = 7
      if (nier8 == 4) info = 8
      if (info /= 0) goto 230
      if (.not.lcv .or. difit > p05*xnorm) goto 20
!
!     Iterative refinement  (if the iteration is converging)
!
      do 210 m = 2,mopt
      fnorm1 = fnorm
      fnorm = 0
      do k = 1,n
         iflag = k
         call sub(n,w(1,2),f,iflag)
         fky = f(k)
         nfcall = nfcall+1
         numf = nfcall/n
         if (iflag < 0) goto 230
         fnorm = max(fnorm,abs(fky))
!
!     Iterative refinement is terminated if it does not give a
!     reduction on residuals
!
         if (fnorm >= fnorm1) then
            fnorm = fnorm1
            goto 20
         end if
         temp = fky/w(k,1)
         do i = 1,n
            w(i,2) = w(i,2)+temp*w(i,k+3)
         end do
      end do
!
!     Compute the norms of the iterate and correction vector
!
      xnorm = 0.d0
      difit = 0.d0
      do i = 1,n
         xnorm = max(xnorm,abs(w(i,2)))
         difit = max(difit,abs(x(i)-w(i,2)))
         x(i) = w(i,2)
      end do
!
!     Stopping criteria for iterative refinement
!
      if (fnorm <= ftol) info = 1
      if (difit <= xtol*xnorm) info = 2
      if (fnorm <= ftol .and. info == 2) info = 3
      if (numf >= maxf .and. info == 0) info = 4
      if (info /= 0) goto 230
  210 continue
      goto 20

  230 if (iflag < 0) info = iflag
      return
      end
