
!     Subroutines to perform the fitting of experimental
!     and theoretical V+A tau spectral function moments
!     with duality violations included
!     Last change: Matthias Jamin, 25.6.2014

!-----------------------------------------------------------------------
      subroutine rtauVpAfit()
      use num_const

      implicit none
      integer :: ierr

      external chisqVpA

      open (unit=7, file='out/minuit_all.out')
      open (unit=8, file='out/fit_all.out')

!     Initialise Minuit
      call mninit(5,6,6)
!     call mninit(5,7,7)
      call mnseti('Fit of alpha_s and non-perturbative parameters')

!     Fit removing D=4 contribution
!     call mnparm( 1,    'astau',  0.31210d0,  0.2d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 2,   'aGGinv',  2.1000d-2,  1.0d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 3, 'rhoD6VpA', -0.56366d-1, 0.4d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 4,  'c8D8VpA', -0.10923d-1, 0.1d-1, 0.d0, 0.d0, ierr)

!     FOPT with D=4
!     call mnparm( 1,    'astau',  0.31777d0,  0.3d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 2,   'aGGinv', -0.21583d-2, 0.2d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 3, 'rhoD6VpA', -0.12953d0,  0.3d-1, 0.d0, 0.d0, ierr)
!     call mnparm( 4,  'c8D8VpA', -0.89766d-1, 0.3d-1, 0.d0, 0.d0, ierr)

!     CIPT with D=4
!     call mnparm( 1,    'astau',  0.33647d0,  0.4d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 2,   'aGGinv', -0.13832d-1, 0.2d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 3, 'rhoD6VpA', -0.77502d-1, 0.2d-1, 0.d0, 0.d0, ierr)
!     call mnparm( 4,  'c8D8VpA', -0.81417d-1, 0.3d-1, 0.d0, 0.d0, ierr)

!     BSPT with D=4
!     call mnparm( 1,    'astau',  0.31757d0,  0.3d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 2,   'aGGinv', -0.13448d-1, 0.2d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 3, 'rhoD6VpA', -0.32649d0,  0.4d-1, 0.d0, 0.d0, ierr)
!     call mnparm( 4,  'c8D8VpA', -0.17804d0,  0.4d-1, 0.d0, 0.d0, ierr)

      call mnparm( 1,    'astau',  0.30978d0,  0.2d-2, 0.d0, 0.d0, ierr)
      call mnparm( 2,   'aGGinv',  2.1000d-2,  1.0d-2, 0.d0, 0.d0, ierr)
      call mnparm( 3, 'rhoD6VpA', -0.12753d0,  0.10d0, 0.d0, 0.d0, ierr)
      call mnparm( 4,  'c8D8VpA',  0.26195d0,  0.30d0, 0.d0, 0.d0, ierr)

      call mnparm( 5,     'c_51',  2.83d2,    1.0d2,  0.d0, 0.d0, ierr)

      call mnparm( 6,     'delV',  3.56d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm( 7,     'gamV',  0.58d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm( 8,     'alpV', -1.92d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm( 9,     'betV',  4.07d0, 0.2d0,  0.d0, 0.d0, ierr)

      call mnparm(10,     'delA',  1.68d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm(11,     'gamA',  1.41d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm(12,     'alpA',  5.16d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm(13,     'betA',  2.13d0, 0.2d0,  0.d0, 0.d0, ierr)


!     Set some fit characteristics
!     call mnexcm(chisqVpA, 'set eps', 1.d-14, 1, ierr, 0)
      call mnexcm(chisqVpA, 'set err', 1.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'set str', 2.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'set pri', 3.d0, 1, ierr, 0)

!     Fix parameters

!     call mnexcm(chisqVpA, 'fix', 1.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 2.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 3.d0, 1, ierr, 0)
!     call mnexcm(chisqVpA, 'fix', 4.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 5.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 6.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 7.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 8.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix', 9.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix',10.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix',11.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix',12.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'fix',13.d0, 1, ierr, 0)
      call mnexcm(chisqVpA, 'call fcn', 0, 0,ierr, 0)

      if (ierr /= 0) then
         write(7,*) "Minuit initialisation error"
         stop
      end if

!     Execute minimisation

      call mnexcm(chisqVpA, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVpA, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVpA, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVpA, 'hesse',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqVpA, 'minos',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqVpA, 'exit',   0.d0, 0, ierr, 0)
!     stop

      call mnintr (chisqVpA,0)

      return
      end subroutine rtauVpAfit
!-----------------------------------------------------------------------

!     Main fit function providing the chi^2

      subroutine chisqVpA(npar,grad,chi2,xval,iflag)
      use params
      use cVA0_const
      use s0s_fit
      use matrix

      implicit none
      integer  :: npar, ierr, iflag
      real(dp) :: chi2, dump, grad(13), xval(13)

      character(len=2) :: scheme = "BS" ! PT-scheme
      integer, parameter :: nevs = Nsum ! # of evs in MILC procedure
      integer, parameter :: Nvar = 3 ! Number of variables in linear
!                                    ! fluctuation analysis

      character(len=8)  :: date
      character(len=10) :: time
      integer  :: k, l, kl, klold, m, n, mn, mnold
      real(dp) :: astau, aGGinv, rhoD6VpA, c8D8VpA, c51,&
                 &deV, kaV, gaV, alV, beV, deA, kaA, gaA, alA, beA
      real(dp) :: momex(Nsum), momth(Nsum)
      real(dp) :: momthVpA(Nsum), dmomthVpA(Nvar,Nsum)
      real(dp) :: covex(Nsum,Nsum), cov0ex(Nsum,Nsum)
      real(dp) :: corex(Nsum,Nsum), cormod(Nsum,Nsum)
      real(dp) :: covinv(Nsum,Nsum), momdif(Nsum)
      real(dp) :: chi2min, chi2n(Nsum), xmin(9)
      real(dp) :: cintDV_VA, cintDVpm_VA, deltaP, sun
      real(dp) :: AmatVpA(Nvar,Nvar), AinvVpA(Nvar,Nvar)
      real(dp) :: BmatVpA(Nvar,Nvar), ADCmatVpA(Nvar,Nsum)
      complex(dc) :: delVpA0FO, delVpA0CI, delVpA0BS
      complex(dc) :: delVpA2FO, delVpA4FO, delVpA68

      common /Minuit/ momex, covex, cov0ex, covinv, chi2min, xmin,&
                     &AmatVpA, AinvVpA, BmatVpA, ADCmatVpA

!     Not to have annoying compiler warnings about 'Unused dummy arguments'
      dump = npar;  dump = grad(1)


!     Initialise some parameters
      if (iflag == 0) then

      call date_and_time(date,time)
      write(8,*) "Starting date/time: ", date, "/", time
      write(8,*) ''

      chi2min = 1.d9

!     Initialise used set of s_0's
      call init_s0s_fit(Nsum,Ns0s,s0s)

!     Initialise moment vectors and covariance matrices (to be safe)
      do k = 1, Nsum
         momex(k) = 0.d0;  momth(k) = 0.d0
         do l = 1, Nsum
            covex(k,l) = 0.d0;  cov0ex(k,l) = 0.d0;  corex(k,l) = 0.d0
         end do
      end do

!     Initialise A-matrices for linear fluctuation analysis
      do m = 1, Nvar
         do n = 1, Nvar
            AmatVpA(m,n) = 0.d0;  AinvVpA(m,n) = 0.d0
            BmatVpA(m,n) = 0.d0
         end do
      end do

      do m = 1, Nvar
         do k = 1, Nsum
            ADCmatVpA(m,k) = 0.d0
         end do
      end do

!     Set up experimental spectral moments (constant during fit)

!     ALEPH MOMENTS
!     Compute experimental moments and covariances for V+A spectral function
      call momex_vpa(momex,covex)

!     Compute experimental correlation matrix
      do kl = 1, Nsum
         do mn = 1, Nsum
             corex(kl,mn) = covex(kl,mn)/&
                           &sqrt(covex(kl,kl))/sqrt(covex(mn,mn))
            cov0ex(kl,mn) = covex(kl,mn)
!     Completely remove correlations between moments
!           if (kl /= mn) then
!              cov0ex(kl,mn) = 0.d0
!           end if
         end do
      end do

!     Remove correlations with R_tau,V+A in Aleph fit
      do kl = 2, Nsum
               cov0ex(1,kl) = 0.d0
               cov0ex(kl,1) = 0.d0
      end do

!     call milc_proc(Nsum,nevs,corex,cormod)

!     Remove cross-correlations between different weights
!     klold = 1;  mnold = 1
!     do k = 1, 3  !  Should be Nmom !!!
!        do kl = klold, klold - 1 + Ns0s(k)
!           do mn = mnold, mnold - 1 + Ns0s(k)
!              cov0ex(kl,mn) = covex(kl,mn)
!              cov0ex(kl,mn+Nsum) = covex(kl,mn+Nsum)
!              cov0ex(kl+Nsum,mn) = covex(kl+Nsum,mn)
!              cov0ex(kl+Nsum,mn+Nsum) = covex(kl+Nsum,mn+Nsum)
!           end do
!        end do
!        klold = kl;  mnold = mn
!     end do

!     Compute inverse of covariance matrices
      call matinv(cov0ex,covinv,Nsum,ierr)
!     write(6,*) "ierr = ", ierr
!     Test inverse
!     call matinvtest(covex,covinv,Nsum)

      end if
!     End of one-time initialisations


!     Set up fit parameters
         astau = xval(1)
        aGGinv = xval(2)
      rhoD6VpA = xval(3)
       c8D8VpA = xval(4)
           c51 = xval(5)
      cV0(5,1,3) = c51
           deV = xval(6)
           gaV = xval(7)
           alV = xval(8)
           beV = xval(9)
           deA = xval(10)
           gaA = xval(11)
           alA = xval(12)
           beA = xval(13)

!     kaV = exp(-deV)
!     kaA = exp(-deA)
      kaV = 0.d0
      kaA = 0.d0

!     Compute theoretical moments and covariances for V spectral function
      call momth_vplusa(Nsum,Ns0s,s0s,scheme, cV0,1.d0,astau,&
                       &aGGinv,rhoD6VpA,c8D8VpA,kaV,gaV,alV,beV,&
                       &kaA,gaA,alA,beA,momthVpA)

!     Difference of experimental and theoretical moments
      do mn = 1, Nsum
         momdif(mn) = momex(mn) - momthVpA(mn)
!        write(*,FMT='(I6,3F12.6)') mn, momex(mn), momthVpA(mn),&
!                                      &momdif(mn)
      end do

!     Compute chi^2

      do k = 1, Nsum
         sun = 0.d0
         do l = 1, Nsum
            sun = sun + momdif(k)*covinv(k,l)*momdif(l)
         end do
         chi2n(k) = sun
!        write(*,FMT='(I6,F16.6)') k, chi2n(k)
      end do

      sun = 0.d0
      do k = 1, Nsum
         sun = sun + chi2n(k)
      end do

      chi2 = sun

      if (chi2 < chi2min) then
         chi2min = chi2
         do k = 1, 15
            xmin(k) = xval(k)
         end do
      end if

      flush(7);  flush(8)

!     Call to user function with iflag = 3
      if (iflag == 3) then

!     Compute linear fluctuation errors for fit parameters

!     Compute derivative of theoretical moments with respect to the
!     theoretical parameters for V spectral function
!     call dmomth_va(Nsum,Ns0s,s0s,scheme, 1,cV0,1.d0,astau,aGGinv,&
!                   &rhoD6V,c8D8V,kaV,gaV,alV,beV,dmomthV)

!     Compute eq. (A6) of paper I
!     do m = 1, Nvar
!        do n = 1, Nvar
!           do k = 1, Nsum
!              do l = 1, Nsum
!                 AmatV(m,n) = AmatV(m,n) +&
!                             &dmomthV(m,k)*covinv(k,l)*dmomthV(n,l)
!              end do
!           end do
!        end do
!     end do

!     call matinv(AmatV,AinvV,Nvar,ierr)

!     Compute eq. (A4) of paper I
!     do m = 1, Nvar
!        do k = 1, Nsum
!           do n = 1, Nvar
!              do l = 1, Nsum
!                 ADCmatV(m,k) = ADCmatV(m,k) +&
!                               &AinvV(m,n)*dmomthV(n,l)*covinv(l,k)
!              end do
!           end do
!        end do
!     end do

!     Compute eq. (A5) of paper I
!     do m = 1, Nvar
!        do n = 1, Nvar
!           do k = 1, Nsum
!              do l = 1, Nsum
!                 BmatV(m,n) = BmatV(m,n) +&
!                             &ADCmatV(m,k)*ADCmatV(n,l)*covex(k,l)
!              end do
!           end do
!        end do
!     end do

!     write(6,*) ''
!     write(6,*) 'Errors with linear fluctuation analysis:'
!     write(6,*) ''
!     write(6,*) "del_alpha_s: ", sqrt(BmatV( 1, 1))
!     write(6,*) "del_C6_V:    ", sqrt(BmatV( 2, 2))
!     write(6,*) "del_C8_V:    ", sqrt(BmatV( 3, 3))

!     write(6,*) ''
!     do m = 1, Nvar
!        write(*,FMT='(7F10.5)') (BmatV(m,n)/&
!                               &sqrt(BmatV(m,m)*BmatV(n,n)), n=1,Nvar)
!     end do

!     Write out fit parameters
      write(8,*) ""
      write(8,*) "alpha_s =   ", astau,    xmin(1)
      write(8,*) "<aGG>_inv = ", aGGinv,   xmin(2)
      write(8,*) "c_6_V+A =   ", rhoD6VpA, xmin(3)
      write(8,*) "c_8_V+A =   ", c8D8VpA,  xmin(4)
      write(8,*) "c_51 =      ", c51,      xmin(5)

      write(8,*) ""
      write(8,*) "FCN = ", chi2, chi2min, iflag

!     call writeVAmomOpal(scheme,astau,aGGinv,rhoD6V,rhoD6A,c8D8V,&
!                        &c8D8A,cV0,kaV,gaV,alV,beV,kaA,gaA,alA,beA)

         write(8,*) ""
         select case (scheme)
         case ("FO")
            write(8,*) "delta_V+A^(0) = ",&
                       &delVpA0FO(1,cV0,stau,stau,astau,5)
         case ("CI")
            write(8,*) "delta_V+A^(0) = ",&
                       &delVpA0CI(1,cV0,stau,stau,astau,5)
         case ("BS")
            write(8,*) "delta_V+A^(0) = ",&
                       &delVpA0BS(1,cV0(5,1,3),stau,astau)
         end select
         write(8,*) "delta_V+A^(2) = ",&
                   &delVpA2FO(1,1,2,stau,stau,astau,2)
         write(8,*) "delta_V+A^(2) = ",&
                   &delVpA2FO(1,1,2,stau,stau,astau,3)
         write(8,*) "delta_V+A^(4) = ",&
                   &delVpA4FO(1,1,2,stau,stau,astau,aGGinv)
         write(8,*) "delta_V+A^(6) = ",&
                   &delVpA68(1,1,2,stau,astau,rhoD6VpA,0.d0)
         write(8,*) "delta_V+A^(8) = ",&
                   &delVpA68(1,1,2,stau,astau,0.d0,c8D8VpA)
         write(8,*) "delta_P+S = ", deltaP(1,stau)
      write(8,*) "deltaDV_V(exact) = ",&
                   &cintDVpm_VA(1,kaV,gaV,alV,beV,stau)/3.d0,&
                   &cintDVpm_VA(1,kaV,gaV,alV,beV,2.89d0)/3.d0
      write(8,*) "deltaDV_V(exact) = ",&
                   &cintDVpm_VA(1,kaV,gaV,alV,beV,2.56D0)/3.d0,&
                   &cintDVpm_VA(1,kaV,gaV,alV,beV,2.25d0)/3.d0
      write(8,*) "deltaDV_A(exact) = ",&
                   &cintDVpm_VA(1,kaA,gaA,alA,beA,stau)/3.d0,&
                   &cintDVpm_VA(1,kaA,gaA,alA,beA,2.89d0)/3.d0
      write(8,*) "deltaDV_A(exact) = ",&
                   &cintDVpm_VA(1,kaA,gaA,alA,beA,2.56D0)/3.d0,&
                   &cintDVpm_VA(1,kaA,gaA,alA,beA,2.25d0)/3.d0
         write(8,*) ""

      call date_and_time(date,time)
      write(8,*) "Ending date/time: ", date, "/", time

      end if

      return
      end subroutine chisqVpA
!-----------------------------------------------------------------------
