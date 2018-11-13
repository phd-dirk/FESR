
!     Main program to fit the vector spectral function to ALEPH data
!     Last change: Matthias Jamin, 12.12.2009

!-----------------------------------------------------------------------

      program fitVV
      use params
      use aleph_vv
      use cVA0_const
      use matrix

      implicit none
      integer, parameter :: N = 140
      integer  :: i, j, Nmin, Nmax, Ns, ierr
      real(dp) :: intsfm2, smin, smax
      real(dp) :: sbin(N), sfm2(N), derr(N), corerr(N,N)
      real(dp) :: correst(N,N), L(N,N)
      real(dp) :: diag(N), rhohat(N)
      real(dp) :: rhoVhat0, sbinj
      complex(dc) :: cintVA0FO, cintVA0CI, cintVA0BS

      external chisqVV

      common/user/ L, sbin, sfm2, derr, smin, smax, rhohat,&
                  &Nmin, Nmax, Ns

!     load ALEPH data
      call aleph_v(sbin,sfm2,derr,corerr)

!     Integrate data to obtain global normalisation factor
      intsfm2 = 0.d0
      do i = 1, N
         intsfm2 = intsfm2 + sfm2(i)*0.025d0
      end do
!     write(*,*) intsfm2  ! intsfm2 = 0.79474684613769564

!     write(6,*) "enter smin and smax "
!     read (5,*) smin, smax
      smin = 1.1d0;  smax = 3.1375d0
      write(*,*) "smin = ",smin, "  smax = ",smax

      Nmin = nint((smin + 0.0125d0)/0.025d0)
      Nmax = nint((smax + 0.0125d0)/0.025d0)
      write(*,*) "Nmin = ", Nmin, "  Nmax = ", Nmax
      Ns = Nmax - Nmin + 1
      write(*,*) "Number of data points = ", Ns

      do i = 1, Ns
         do j = 1, Ns
            correst(i,j) = corerr(i+Nmin-1,j+Nmin-1)
         end do
      end do

!     Cholesky decomposition of corerr

      call choldc(correst,Ns,N,diag)
      do i = 1, Ns
         do j = 1, i-1
            L(i,j) = correst(i,j)
         end do
         L(i,i) = diag(i)
      end do

!     Prepare radiative correction for spectral function
      open (unit=9, file='rhohat.dat')
      do i = 1, Ns
         j = i + Nmin - 1
         sbinj = sbin(j)  ! Why is this necessary?
!     Choose between FOPT, CIPT or Borel sum   
!        rhohat(j) = rhoVhat0(3,sbinj,sbinj,astauBJ,5)     ! FOPT explicitly
!        rhohat(j) = real(cintVA0FO(0,cV0,sbinj,sbinj,astauBJ,5),dp) ! FOPT
!        rhohat(j) = real(cintVA0CI(0,cV0,sbinj,sbinj,astauBJ,5),dp) ! CIPT
         rhohat(j) = real(cintVA0BS(0,cV0(5,1,3),sbinj,astauBJ),dp)  ! Borel sum
         write(9,*) sbinj, rhohat(j)
      end do
      close (9)

!     MINUIT

!     Initialise Minuit                                                     
!                                                                       
      call mninit(5,6,8)
      call mnseti('Fit of the ALEPH vector spectral function')
!     call mnexcm(chisqVV, 'set print', 1, 1,ierr, 0)
      call mnparm(1, 'astau',  0.3156d0, 0.002d0, 0.d0, 0.d0, ierr)
      call mnparm(2, 'kappa', -0.0172d0, 0.005d0, 0.d0, 0.d0, ierr)
      call mnparm(3, 'gamma',  0.149d0 , 0.05d0 , 0.d0, 5.d0, ierr)
      call mnparm(4, 'alpha',  4.24d0  , 0.5d0  , 0.d0, 0.d0, ierr)
      call mnparm(5, 'beta' , -2.08d0  , 0.5d0  , 0.d0, 0.d0, ierr)
      call mnparm(6, 'a4',     0.d0    , 0.2d0  , 0.d0, 0.d0, ierr)
      if (ierr /= 0) then
         write(*,*) "Minuit initialisation error"
         stop
      end if
!     call mnexcm(chisqVV, 'set eps', 1.d-14, 1, ierr, 0)
      call mnexcm(chisqVV, 'set str', 2.d0, 1, ierr, 0)
      call mnexcm(chisqVV, 'set err', 1.d0, 1, ierr, 0)

!     Fix parameters

      call mnexcm(chisqVV, 'fix', 1.d0, 1, ierr, 0)
      call mnexcm(chisqVV, 'fix', 6.d0, 1, ierr, 0)
      call mnexcm(chisqVV, 'call fcn', 0, 0,ierr, 0)

!     Execute minimisation

      call mnexcm(chisqVV, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVV, 'hesse',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqVV, 'minos',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqVV, 'stop',   0.d0, 0, ierr, 0)

!     stop 

      call mnintr (chisqVV,0)

      end

!-----------------------------------------------------------------------

!     Main fit function providing the chi^2

      subroutine chisqVV(npar,grad,chi2,xval,iflag)
      use params
      use matrix

      implicit none
      integer  :: npar, iflag
      real(dp) :: chi2, grad(6), xval(6)

      integer,  parameter :: M = 140
      real(dp), parameter :: Norm = 1.2195850385647344d0
!     Norm = M_tau^2/(6|V_ud|^2 S_EW) R_tau,V/intsfm2 =
!            1.77684^2/(6 0.97425^2 1.0198) 1.783/0.79474684613769564
!     !!! Compute error for Norm !!!

      integer     :: i , j, Nmin, Nmax, Ns
      real(dp)    :: smin, smax, sun
      real(dp)    :: sbin(M), sfm2(M), derr(M)
      real(dp)    :: L(M,M), v(M), rhohat(M)
      real(dp)    :: rhoVhat0, sbinj
!     complex(dc) :: cintVA0BS

      common/user/ L, sbin, sfm2, derr, smin, smax, rhohat,&
                  &Nmin, Nmax, Ns

!     L contains the result of the Cholesky decomposition of the correlation
!     matrix (lower trangular part). Now get L^{-1}v, with v the vector of
!     differences between the theoretical and experimental spectral function
!     (normalized as Aleph's v1), multiplied by the (s-dependent) normalisation
!     factor (1-sbin(i)/m_tau^2)^2(1+2sbin(i)/m_tau^2)/Norm, weighted by
!     1/derr(i) (to convert from inverse correlation to inverse covariance
!     matrix in chi2), and by a factor 10 because corerr is given in percents
!     (see aleph_v):                    

      do i = 1, Ns
         j = i + Nmin - 1
         sbinj = sbin(j)
!        rhohat(j) = rhoVhat0(3,sbinj,sbinj,xval(1),5)   ! FOPT
!     Presently too slow to be used in the fit
!        rhohat(j) = real(cintVA0BS(0,cV0(5,1,3),sbinj,xval(1)),dp) ! Borel sum
         v(i) = 10.d0*( ( 0.5d0*( 1.d0 + rhohat(j) ) +&
               &xval(6)/sbinj**2 +&
               &2.d0*pi**2*xval(2)*exp(-xval(3)*sbinj) *&
               &sin(xval(4)+xval(5)*sbinj) ) *&
               &(1.d0-sbinj/stau)**2*(1.d0+2.d0*sbinj/stau)/Norm-&
               &sfm2(j))/derr(j)
      end do 
      call lbksb(L,Ns,M,v) 

!     Compute chi^2

      sun = 0.d0
      do i = 1, Ns
         sun = sun + v(i)*v(i)
      end do

      chi2 = sun

!     Call to user function with iflag = 3

      if (iflag == 3) then
         write(*,*) " alpha_s = ", xval(1), " kappa = ", xval(2),&
                   &" gamma = ", xval(3), " alpha = ",  xval(4),&
                   &" beta = ", xval(5), " a4 = ", xval(6)
      end if
      return
      end
!-----------------------------------------------------------------------
