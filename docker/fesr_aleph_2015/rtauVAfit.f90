
!     Subroutines to perform the fitting of experimental
!     and theoretical tau spectral function moments
!     with duality violations included
!     Last change: Matthias Jamin, 26.2.2012

!-----------------------------------------------------------------------
      subroutine rtauVAfit()
      use num_const

      implicit none
      integer :: ierr

      external chisqVA

      open (unit=7, file='out/minuit_all.out')
      open (unit=8, file='out/fit_all.out')

!     Initialise Minuit
      call mninit(5,6,6)
!     call mninit(5,7,7)
      call mnseti('Fit of alpha_s and non-perturbative parameters')

      call mnparm( 1,  'astau', 0.31113d0, 0.5d-2, 0.d0, 0.d0, ierr)
!     call mnparm( 2, 'aGGinv', 0.0000d-2, 1.0d-2, 0.d0, 0.d0, ierr)
      call mnparm( 2, 'aGGinv', 2.1000d-2, 1.0d-2, 0.d0, 0.d0, ierr)
      call mnparm( 3, 'rhoD6V', 0.53739d0, 0.10d0, 0.d0, 0.d0, ierr)
      call mnparm( 4, 'rhoD6A',-0.47785d0, 0.10d0, 0.d0, 0.d0, ierr)
      call mnparm( 5,  'c8D8V', 0.93364d0, 0.10d0, 0.d0, 0.d0, ierr)
      call mnparm( 6,  'c8D8A',-0.32840d0, 0.10d0, 0.d0, 0.d0, ierr)
      call mnparm( 7,   'c_51', 2.83d2,    1.0d2,  0.d0, 0.d0, ierr)

      call mnparm( 8,   'delV', 3.6612d0,  0.3d0,  0.d0, 0.d0, ierr)
      call mnparm( 9,   'gamV', 0.42755d0, 0.2d0,  0.d0, 0.d0, ierr)
      call mnparm(10,   'alpV',-1.2544d0,  0.2d0,  0.d0, 0.d0, ierr)
      call mnparm(11,   'betV', 3.7583d0,  0.1d0,  0.d0, 0.d0, ierr)

      call mnparm(12,   'delA', 2.1690d0,  0.1d0,  0.d0, 0.d0, ierr)
      call mnparm(13,   'gamA', 1.2657d0,  0.1d0,  0.d0, 0.d0, ierr)
      call mnparm(14,   'alpA',-2.7239d0,  0.2d0,  0.d0, 0.d0, ierr)
      call mnparm(15,   'betA', 3.0696d0,  0.1d0,  0.d0, 0.d0, ierr)


!     Set some fit characteristics
!     call mnexcm(chisqVA, 'set eps', 1.d-14, 1, ierr, 0)
      call mnexcm(chisqVA, 'set err', 1.d0, 1, ierr, 0)
      call mnexcm(chisqVA, 'set str', 2.d0, 1, ierr, 0)
      call mnexcm(chisqVA, 'set pri', 3.d0, 1, ierr, 0)

!     Fix parameters

!     call mnexcm(chisqVA, 'fix', 1.d0, 1, ierr, 0)
      call mnexcm(chisqVA, 'fix', 2.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix', 3.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix', 4.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix', 5.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix', 6.d0, 1, ierr, 0)
      call mnexcm(chisqVA, 'fix', 7.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix', 8.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix', 9.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix',10.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix',11.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix',12.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix',13.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix',14.d0, 1, ierr, 0)
!     call mnexcm(chisqVA, 'fix',15.d0, 1, ierr, 0)
      call mnexcm(chisqVA, 'call fcn', 0, 0,ierr, 0)

      if (ierr /= 0) then
         write(7,*) "Minuit initialisation error"
         stop
      end if

!     Execute minimisation

!     call mnexcm(chisqVA, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVA, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVA, 'migrad', 0.d0, 0, ierr, 0)
!     call mnexcm(chisqVA, 'hesse',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqVA, 'minos',  0.d0, 0, ierr, 0)
!     call mnexcm(chisqVA, 'exit',   0.d0, 0, ierr, 0)
!     stop

      call mnintr (chisqVA,0)

      return
      end subroutine rtauVAfit
!-----------------------------------------------------------------------

!     Main fit function providing the chi^2

      subroutine chisqVA(npar,grad,chi2,xval,iflag)
      use params
      use cVA0_const
      use s0s_fit
      use matrix

      implicit none
      integer  :: npar, ierr, iflag
      real(dp) :: chi2, dump, grad(15), xval(15)

      character(len=2) :: scheme = "FO" ! PT-scheme
      integer, parameter :: nevs = 2*Nsum ! # of evs in MILC procedure
      integer, parameter :: Nvar = 13 ! Number of variables in linear
!                                     ! fluctuation analysis

      character(len=8)  :: date
      character(len=10) :: time
      integer  :: k, l, kl, klold, m, n, mn, mnold
      real(dp) :: astau, aGGinv, rhoD6V, rhoD6A, c8D8V, c8D8A, c51,&
                 &deV, kaV, gaV, alV, beV, deA, kaA, gaA, alA, beA
      real(dp) :: momex(2*Nsum), momth(2*Nsum)
      real(dp) :: dmomthV(7,Nsum), dmomthA(7,Nsum)
      real(dp) :: momthV(Nsum), momthA(Nsum), dmomthVA(Nvar,2*Nsum)
      real(dp) :: covex(2*Nsum,2*Nsum), cov0ex(2*Nsum,2*Nsum)
      real(dp) :: corex(2*Nsum,2*Nsum), cormod(2*Nsum,2*Nsum)
      real(dp) :: covinv(2*Nsum,2*Nsum), momdif(2*Nsum)
      real(dp) :: chi2min, chi2n(2*Nsum), xmin(15)
      real(dp) :: cintDV_VA, cintDVpm_VA, sun
      real(dp) :: AmatVA(Nvar,Nvar), AinvVA(Nvar,Nvar)
      real(dp) :: BmatVA(Nvar,Nvar), ADCmatVA(Nvar,2*Nsum)
      complex(dc) :: delVA0FO, delVA0CI, delVA0BS
      complex(dc) :: delVA2FO, delVA4FO, delVA68, deltaP

      common /Minuit/ momex, covex, cov0ex, covinv, chi2min, xmin,&
                     &AmatVA, AinvVA, BmatVA, ADCmatVA

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
            AmatVA(m,n) = 0.d0;  AinvVA(m,n) = 0.d0;  BmatVA(m,n) = 0.d0
         end do
      end do

      do m = 1, Nvar
         do k = 1, 2*Nsum
            ADCmatVA(m,k) = 0.d0
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

      do kl = 1, 2*Nsum
         write(*,*) "kl = ", kl, "   ", momex(kl), "   ",&
                    &dsqrt(covex(kl,kl))
      end do
!     write(*,*) ""
!     write(*,*) covex(1,55)/dsqrt(covex(1,1)*covex(55,55))
      
!     Compute experimental correlation matrix
      do kl = 1, 2*Nsum
         do mn = 1, 2*Nsum
            corex(kl,mn) = covex(kl,mn)/&
                          &sqrt(covex(kl,kl))/sqrt(covex(mn,mn))
!           cov0ex(kl,mn) = covex(kl,mn)
!     Completely remove correlations between moments
!           if (kl /= mn) then
!              covex(kl,mn) = 0.d0
!           end if
         end do
      end do

      call milc_proc(2*Nsum,nevs,corex,cormod)

!     Remove cross-correlations between different weights
      klold = 1;  mnold = 1
      do k = 1, 3  !  Should be Nmom !!!
         do kl = klold, klold - 1 + Ns0s(k)
            do mn = mnold, mnold - 1 + Ns0s(k)
               cov0ex(kl,mn) = covex(kl,mn)
               cov0ex(kl,mn+Nsum) = covex(kl,mn+Nsum)
               cov0ex(kl+Nsum,mn) = covex(kl+Nsum,mn)
               cov0ex(kl+Nsum,mn+Nsum) = covex(kl+Nsum,mn+Nsum)
            end do
         end do
         klold = kl;  mnold = mn
      end do

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
      rhoD6A = xval(4)
       c8D8V = xval(5)
       c8D8A = xval(6)
         c51 = xval(7)
      cV0(5,1,3) = c51
         deV = xval(8)
         gaV = xval(9)
         alV = xval(10)
         beV = xval(11)
         deA = xval(12)
         gaA = xval(13)
         alA = xval(14)
         beA = xval(15)

      kaV = exp(-deV)
      kaA = exp(-deA)


!     Compute theoretical moments and covariances for V spectral function
      call momth_va(Nsum,Ns0s,s0s,scheme, 1,cV0,1.d0,astau,aGGinv,&
                   &rhoD6V,c8D8V,kaV,gaV,alV,beV,momthV)

!     Compute theoretical moments and covariances for A spectral function
      call momth_va(Nsum,Ns0s,s0s,scheme,-1,cV0,1.d0,astau,aGGinv,&
                   &rhoD6A,c8D8A,kaA,gaA,alA,beA,momthA)


!     Difference of experimental and theoretical moments
      do mn = 1, Nsum
         momdif(mn)      = momex(mn)      - momthV(mn)
         momdif(mn+Nsum) = momex(mn+Nsum) - momthA(mn)
!        write(*,FMT='(I6,3F12.6)') mn, momex(mn), momthV(mn),&
!                                      &momdif(mn)
!        write(*,FMT='(I6,3F12.6)') mn, momex(mn+Nsum), momthA(mn),&
!                                      &momdif(mn+Nsum)
!        write(*,FMT='(I6,2F12.6)') mn, momex(mn)+momex(mn+Nsum),&
!                                      &momthV(mn)+momthA(mn)
      end do

!     Compute chi^2

      do k = 1, 2*Nsum
         sun = 0.d0
         do l = 1, 2*Nsum
            sun = sun + momdif(k)*covinv(k,l)*momdif(l)
         end do
         chi2n(k) = sun
!        write(*,FMT='(I6,F16.6)') k, chi2n(k)
      end do

      sun = 0.d0
      do k = 1, 2*Nsum
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

!     Compute derivative of theoretical moments with respect to the
!     theoretical parameters for A spectral function
      call dmomth_va(Nsum,Ns0s,s0s,scheme,-1,cV0,1.d0,astau,aGGinv,&
                    &rhoD6A,c8D8A,kaA,gaA,alA,beA,dmomthA)

!     set up dmomthVA
      do k = 1, Nsum
         dmomthVA(1,k     ) = dmomthV(1,k)
         dmomthVA(1,k+Nsum) = dmomthA(1,k)
      end do

      do m = 2, 7
         do k = 1, Nsum
            dmomthVA(m,k     ) = dmomthV(m,k)
            dmomthVA(m,k+Nsum) = 0.d0
         end do
      end do
         
      do m = 2, 7
         do k = 1, Nsum
            dmomthVA(m+6,k     ) = 0.d0
            dmomthVA(m+6,k+Nsum) = dmomthA(m,k)
         end do
      end do

!     Compute eq. (A6) of paper I
      do m = 1, Nvar
         do n = 1, Nvar
            do k = 1, 2*Nsum
               do l = 1, 2*Nsum
                  AmatVA(m,n) = AmatVA(m,n) +&
                               &dmomthVA(m,k)*covinv(k,l)*dmomthVA(n,l)
               end do
            end do
         end do
      end do

      call matinv(AmatVA,AinvVA,Nvar,ierr)

!     Compute eq. (A4) of paper I
      do m = 1, Nvar
         do k = 1, 2*Nsum
            do n = 1, Nvar
               do l = 1, 2*Nsum
                  ADCmatVA(m,k) = ADCmatVA(m,k) +&
                                 &AinvVA(m,n)*dmomthVA(n,l)*covinv(l,k)
               end do
            end do
         end do
      end do

!     Compute eq. (A5) of paper I
      do m = 1, Nvar
         do n = 1, Nvar
            do k = 1, 2*Nsum
               do l = 1, 2*Nsum
                  BmatVA(m,n) = BmatVA(m,n) +&
                               &ADCmatVA(m,k)*ADCmatVA(n,l)*covex(k,l)
               end do
            end do
         end do
      end do

      write(6,*) ''
      write(6,*) 'Errors with linear fluctuation analysis:'
      write(6,*) ''
      write(6,*) "del_alpha_s: ", sqrt(BmatVA( 1, 1))
      write(6,*) "del_C6_V:    ", sqrt(BmatVA( 2, 2))
      write(6,*) "del_C8_V:    ", sqrt(BmatVA( 3, 3))
      write(6,*) "del_delta_V: ", sqrt(BmatVA( 4, 4))
      write(6,*) "del_gamma_V: ", sqrt(BmatVA( 5, 5))
      write(6,*) "del_alpha_V: ", sqrt(BmatVA( 6, 6))
      write(6,*) "del_beta_V:  ", sqrt(BmatVA( 7, 7))
      write(6,*) "del_C6_A:    ", sqrt(BmatVA( 8, 8))
      write(6,*) "del_C8_A:    ", sqrt(BmatVA( 9, 9))
      write(6,*) "del_delta_A: ", sqrt(BmatVA(10,10))
      write(6,*) "del_gamma_A: ", sqrt(BmatVA(11,11))
      write(6,*) "del_alpha_A: ", sqrt(BmatVA(12,12))
      write(6,*) "del_beta_A:  ", sqrt(BmatVA(13,13))

      write(6,*) ''
      do m = 1, Nvar
         write(*,FMT='(13F10.5)') (BmatVA(m,n)/&
                             &sqrt(BmatVA(m,m)*BmatVA(n,n)), n=1,Nvar)
      end do

!     Write out fit parameters
      write(8,*) ""
      write(8,*) "alpha_s =   ", astau,  xmin(1)
      write(8,*) "<aGG>_inv = ", aGGinv, xmin(2)
      write(8,*) "rho_V =     ", rhoD6V, xmin(3)
      write(8,*) "rho_A =     ", rhoD6A, xmin(4)
      write(8,*) "c_8,V =     ", c8D8V,  xmin(5)
      write(8,*) "c_8,A =     ", c8D8A,  xmin(6)
      write(8,*) "c_51 =      ", c51,    xmin(7)
      write(8,*) "delta_V =   ", deV,    xmin(8)
      write(8,*) "gamma_V =   ", gaV,    xmin(9)
      write(8,*) "alpha_V =   ", alV,    xmin(10)
      write(8,*) "beta_V =    ", beV,    xmin(11)
      write(8,*) "delta_A =   ", deA,    xmin(12)
      write(8,*) "gamma_A =   ", gaA,    xmin(13)
      write(8,*) "alpha_A =   ", alA,    xmin(14)
      write(8,*) "beta_A =    ", beA,    xmin(15)

      write(8,*) ""
      write(8,*) "chi^2 = ", chi2, chi2min, iflag

!     call writeVAmomOpal(scheme,astau,aGGinv,rhoD6V,rhoD6A,c8D8V,&
!                        &c8D8A,cV0,kaV,gaV,alV,beV,kaA,gaA,alA,beA)

!     Call to user function with iflag = 3
!     if (iflag == 3) then

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
         write(8,*) "delta_A^(4) = ",&
                   &delVA4FO(1,-1,1,2,stau,stau,astau,aGGinv)
         write(8,*) "delta_V^(6) = ",&
                   &delVA68(3, 1,1,2,stau,astau,rhoD6V,0.d0)
         write(8,*) "delta_A^(6) = ",&
                   &delVA68(3,-1,1,2,stau,astau,rhoD6A,0.d0)
         write(8,*) "delta_V^(8) = ",&
                   &delVA68(3, 1,1,2,stau,astau,0.d0,c8D8V)
         write(8,*) "delta_A^(8) = ",&
                   &delVA68(3,-1,1,2,stau,astau,0.d0,c8D8A)
         write(8,*) "delta_P+S = ", deltaP(1,stau)
      write(8,*) "DV_V(exact) = ",&
                   &cintDVpm_VA(1,kaV,gaV,alV,beV,stau),&
                   &cintDVpm_VA(3,kaV,gaV,alV,beV,stau)
      write(8,*) "DV_V(integ) = ",&
                   &cintDV_VA(1,kaV,gaV,alV,beV,stau),&
                   &cintDV_VA(3,kaV,gaV,alV,beV,stau)
      write(8,*) "DV_A(exact) = ",&
                   &cintDVpm_VA(1,kaA,gaA,alA,beA,stau),&
                   &cintDVpm_VA(3,kaA,gaA,alA,beA,stau)
      write(8,*) "DV_A(integ) = ",&
                   &cintDV_VA(1,kaA,gaA,alA,beA,stau),&
                   &cintDV_VA(3,kaA,gaA,alA,beA,stau)
         write(8,*) ""

      call date_and_time(date,time)
      write(8,*) "Ending date/time: ", date, "/", time

      end if

      return
      end subroutine chisqVA
!-----------------------------------------------------------------------
