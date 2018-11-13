
!     Subroutines to perform the fitting of experimental
!     and theoretical tau spectral function moments
!     with duality violations included
!     Last change: Matthias Jamin, 26.2.2012

!-----------------------------------------------------------------------
      subroutine rtauVfit()
      use num_const

      implicit none
      integer :: ierr

      external chisqV

      open (unit=7, file='out/minuit_all.out')
      open (unit=8, file='out/fit_all.out')

!     Initialise Minuit
      call mninit(5,6,6)
!     call mninit(5,7,7)
      call mnseti('Fit of alpha_s and non-perturbative parameters')

      call mnparm( 1,  'astau', 0.32341d0, 1.6d-2, 0.d0, 0.d0, ierr)
      call mnparm( 2, 'aGGinv', 2.1000d-2, 1.0d-2, 0.d0, 0.d0, ierr)
      call mnparm( 3, 'rhoD6V', 0.00000d0, 0.10d0, 0.d0, 0.d0, ierr)
      call mnparm( 4,  'c8D8V', 0.00000d0, 0.20d0, 0.d0, 0.d0, ierr)
      call mnparm( 5,   'c_51', 2.83d2,    1.0d2,  0.d0, 0.d0, ierr)

      call mnparm( 6,   'delV', 4.2096d0,  0.6d0,  0.d0, 0.d0, ierr)
      call mnparm( 7,   'gamV', 0.11703d0, 0.4d0,  0.d0, 0.d0, ierr)
      call mnparm( 8,   'alpV',-0.48473d0, 0.8d0,  0.d0, 0.d0, ierr)
      call mnparm( 9,   'betV', 3.3786d0,  0.4d0,  0.d0, 0.d0, ierr)


!     Set some fit characteristics
!     call mnexcm(chisqV, 'set eps', 1.d-14, 1, ierr, 0)
      call mnexcm(chisqV, 'set err', 1.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'set str', 2.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'set pri', 3.d0, 1, ierr, 0)

!     Fix parameters

!     call mnexcm(chisqV, 'fix', 1.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'fix', 2.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'fix', 3.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'fix', 4.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'fix', 5.d0, 1, ierr, 0)
!     call mnexcm(chisqV, 'fix', 6.d0, 1, ierr, 0)
!     call mnexcm(chisqV, 'fix', 7.d0, 1, ierr, 0)
!     call mnexcm(chisqV, 'fix', 8.d0, 1, ierr, 0)
!     call mnexcm(chisqV, 'fix', 9.d0, 1, ierr, 0)
      call mnexcm(chisqV, 'call fcn', 0, 0,ierr, 0)

      if (ierr /= 0) then
         write(7,*) "Minuit initialisation error"
         stop
      end if

!     Execute minimisation

      call mnexcm(chisqV, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqV, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqV, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqV, 'hesse',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqV, 'minos',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqV, 'exit',   0.d0, 0, ierr, 0)
!     stop

      call mnintr (chisqV,0)

      return
      end subroutine rtauVfit
!-----------------------------------------------------------------------

!     Main fit function providing the chi^2

      subroutine chisqV(npar,grad,chi2,xval,iflag)
      use params
      use cVA0_const
      use s0s_fit
      use matrix

      implicit none
      integer  :: npar, ierr, iflag
      real(dp) :: chi2, dump, grad(15), xval(15)

      character(len=2) :: scheme = "FO" ! PT-scheme
      integer, parameter :: nevs = 2*Nsum ! # of evs in MILC procedure
      integer, parameter :: Nvar = 7 ! Number of variables in linear
!                                    ! fluctuation analysis

      character(len=8)  :: date
      character(len=10) :: time
      integer  :: k, l, kl, klold, m, n, mn, mnold
      real(dp) :: astau, aGGinv, rhoD6V, c8D8V, c51,&
                 &deV, kaV, gaV, alV, beV
      real(dp) :: momex(2*Nsum), momth(2*Nsum)
      real(dp) :: momthV(Nsum), dmomthV(Nvar,Nsum)
      real(dp) :: covex(2*Nsum,2*Nsum), cov0ex(2*Nsum,2*Nsum)
      real(dp) :: corex(2*Nsum,2*Nsum), cormod(2*Nsum,2*Nsum)
      real(dp) :: covinv(2*Nsum,2*Nsum), momdif(2*Nsum)
      real(dp) :: chi2min, chi2n(2*Nsum), xmin(15)
      real(dp) :: cintDV_VA, cintDVpm_VA, sun
      real(dp) :: AmatV(Nvar,Nvar), AinvV(Nvar,Nvar)
      real(dp) :: BmatV(Nvar,Nvar), ADCmatV(Nvar,Nsum)
      complex(dc) :: delVA0FO, delVA0CI, delVA0BS
      complex(dc) :: delVA2FO, delVA4FO, delVA68, deltaP

      common /Minuit/ momex, covex, cov0ex, covinv, chi2min, xmin,&
                     &AmatV, AinvV, BmatV, ADCmatV

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
      do k = 1, 2*Nsum
         momex(k) = 0.d0;  momth(k) = 0.d0
         do l = 1, 2*Nsum
            covex(k,l) = 0.d0;  cov0ex(k,l) = 0.d0;  corex(k,l) = 0.d0
         end do
      end do

!     Initialise A-matrices for linear fluctuation analysis
      do m = 1, Nvar
         do n = 1, Nvar
            AmatV(m,n) = 0.d0;  AinvV(m,n) = 0.d0;  BmatV(m,n) = 0.d0
         end do
      end do

      do m = 1, Nvar
         do k = 1, Nsum
            ADCmatV(m,k) = 0.d0
         end do
      end do

!     Set up experimental spectral moments (constant during fit)

!     ALEPH MOMENTS ! (Subroutine needs to be updated)
!     Compute experimental moments and covariances for V spectral function
!     call    momex_v(momexV,covexV)
!     call relmomex_v(momexV,covexV)

!     Compute experimental moments and covariances for A spectral function
!     call    momex_a(momexA,covexA)
!     call relmomex_a(momexA,covexA)

!     OPAL MOMENTS
!     Compute experimental moments and covariances for V/A spectral function
      call momex_va(Nsum,Ns0s,s0s,momex,covex)

!     Compute experimental correlation matrix
      do kl = 1, 2*Nsum
         do mn = 1, 2*Nsum
             corex(kl,mn) = covex(kl,mn)/&
                           &sqrt(covex(kl,kl))/sqrt(covex(mn,mn))
            cov0ex(kl,mn) = covex(kl,mn)
!     Completely remove correlations between moments
!           if (kl /= mn) then
!              cov0ex(kl,mn) = 0.d0
!           end if
         end do
      end do

      call milc_proc(2*Nsum,nevs,corex,cormod)

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
      call matinv(cov0ex,covinv,2*Nsum,ierr)
!     write(8,*) "ierr = ", ierr
!     Test inverse
!     call matinvtest(covex,covinv,2*Nsum)

      end if
!     End of one-time initialisations


!     Set up fit parameters
       astau = xval(1)
      aGGinv = xval(2)
      rhoD6V = xval(3)
       c8D8V = xval(4)
         c51 = xval(5)
      cV0(5,1,3) = c51
         deV = xval(6)
         gaV = xval(7)
         alV = xval(8)
         beV = xval(9)

      kaV = exp(-deV)


!     Compute theoretical moments and covariances for V spectral function
      call momth_va(Nsum,Ns0s,s0s,scheme, 1,cV0,1.d0,astau,aGGinv,&
                   &rhoD6V,c8D8V,kaV,gaV,alV,beV,momthV)

!     Difference of experimental and theoretical moments
      do mn = 1, Nsum
         momdif(mn) = momex(mn) - momthV(mn)
!        write(*,FMT='(I6,3F12.6)') mn, momex(mn), momthV(mn),&
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
      call dmomth_va(Nsum,Ns0s,s0s,scheme, 1,cV0,1.d0,astau,aGGinv,&
                    &rhoD6V,c8D8V,kaV,gaV,alV,beV,dmomthV)

!     Compute eq. (A6) of paper I
      do m = 1, Nvar
         do n = 1, Nvar
            do k = 1, Nsum
               do l = 1, Nsum
                  AmatV(m,n) = AmatV(m,n) +&
                              &dmomthV(m,k)*covinv(k,l)*dmomthV(n,l)
               end do
            end do
         end do
      end do

      call matinv(AmatV,AinvV,Nvar,ierr)

!     Compute eq. (A4) of paper I
      do m = 1, Nvar
         do k = 1, Nsum
            do n = 1, Nvar
               do l = 1, Nsum
                  ADCmatV(m,k) = ADCmatV(m,k) +&
                                &AinvV(m,n)*dmomthV(n,l)*covinv(l,k)
               end do
            end do
         end do
      end do

!     Compute eq. (A5) of paper I
      do m = 1, Nvar
         do n = 1, Nvar
            do k = 1, Nsum
               do l = 1, Nsum
                  BmatV(m,n) = BmatV(m,n) +&
                              &ADCmatV(m,k)*ADCmatV(n,l)*covex(k,l)
               end do
            end do
         end do
      end do

      write(6,*) ''
      write(6,*) 'Errors with linear fluctuation analysis:'
      write(6,*) ''
      write(6,*) "del_alpha_s: ", sqrt(BmatV( 1, 1))
      write(6,*) "del_C6_V:    ", sqrt(BmatV( 2, 2))
      write(6,*) "del_C8_V:    ", sqrt(BmatV( 3, 3))
      write(6,*) "del_delta_V: ", sqrt(BmatV( 4, 4))
      write(6,*) "del_gamma_V: ", sqrt(BmatV( 5, 5))
      write(6,*) "del_alpha_V: ", sqrt(BmatV( 6, 6))
      write(6,*) "del_beta_V:  ", sqrt(BmatV( 7, 7))

      write(6,*) ''
      do m = 1, Nvar
         write(*,FMT='(7F10.5)') (BmatV(m,n)/&
                                &sqrt(BmatV(m,m)*BmatV(n,n)), n=1,Nvar)
      end do

!     Write out fit parameters
      write(8,*) ""
      write(8,*) "alpha_s =   ", astau,  xmin(1)
      write(8,*) "<aGG>_inv = ", aGGinv, xmin(2)
      write(8,*) "rho_V =     ", rhoD6V, xmin(3)
      write(8,*) "c_8,V =     ", c8D8V,  xmin(4)
      write(8,*) "c_51 =      ", c51,    xmin(5)
      write(8,*) "delta_V =   ", deV,    xmin(6)
      write(8,*) "gamma_V =   ", gaV,    xmin(7)
      write(8,*) "alpha_V =   ", alV,    xmin(8)
      write(8,*) "beta_V =    ", beV,    xmin(9)

      write(8,*) ""
      write(8,*) "chi^2 = ", chi2, chi2min, iflag

!     call writeVAmomOpal(scheme,astau,aGGinv,rhoD6V,rhoD6A,c8D8V,&
!                        &c8D8A,cV0,kaV,gaV,alV,beV,kaA,gaA,alA,beA)

         write(8,*) ""
         select case (scheme)
         case ("FO")
            write(8,*) "delta_V/A^(0) = ",&
                       &delVA0FO(1,cV0,stau,stau,astau,5)
         case ("CI")
            write(8,*) "delta_V/A^(0) = ",&
                       &delVA0CI(1,cV0,stau,stau,astau,5)
         case ("BS")
            write(8,*) "delta_V/A^(0) = ",&
                       &delVA0BS(1,cV0(5,1,3),stau,astau)
         end select
         write(8,*) "delta_V^(2) = ",&
                   &delVA2FO(1, 1,1,2,stau,stau,astau,2)
         write(8,*) "delta_V^(2) = ",&
                   &delVA2FO(1, 1,1,2,stau,stau,astau,3)
         write(8,*) "delta_A^(2) = ",&
                   &delVA2FO(1,-1,1,2,stau,stau,astau,2)
         write(8,*) "delta_A^(2) = ",&
                   &delVA2FO(1,-1,1,2,stau,stau,astau,3)
         write(8,*) "delta_V^(4) = ",&
                   &delVA4FO(1, 1,1,2,stau,stau,astau,aGGinv)
         write(8,*) "delta_V^(6) = ",&
                   &delVA68(3, 1,1,2,stau,astau,rhoD6V,0.d0)
         write(8,*) "delta_V^(8) = ",&
                   &delVA68(3, 1,1,2,stau,astau,0.d0,c8D8V)
         write(8,*) "delta_P+S = ", deltaP(1,stau)
      write(8,*) "DV_V(exact) = ",&
                   &cintDVpm_VA(1,kaV,gaV,alV,beV,stau),&
                   &cintDVpm_VA(2,kaV,gaV,alV,beV,stau)
      write(8,*) "DV_V(integ) = ",&
                   &cintDV_VA(1,kaV,gaV,alV,beV,stau),&
                   &cintDV_VA(2,kaV,gaV,alV,beV,stau)
         write(8,*) ""

      call date_and_time(date,time)
      write(8,*) "Ending date/time: ", date, "/", time

      end if

      return
      end subroutine chisqV
!-----------------------------------------------------------------------
