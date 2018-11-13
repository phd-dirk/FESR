
!     Compute experimental moments with correlations from OPAL data
!     Employing 2011 update of the OPAL data by Kim Maltman
!     Last change: Matthias Jamin, 20.2.2012

!-----------------------------------------------------------------------

!     Moment integrals of V/A spectral function

      subroutine momex_va(Nsum,Ns0s,s0s,mom,covmom)
      use params
!     use opal98
      use weights_rho
      use weights_used

      implicit none
      integer,  intent(in) :: Nsum ! Total number of s_0's
      integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
      real(dp), dimension(Nmom,0:140) :: s0s ! Sets of used s_0's
      real(dp), intent(out) :: mom(2*Nsum), covmom(2*Nsum,2*Nsum)

      integer, parameter :: Ndat = 98 ! Number of new Opal data points
      integer  :: Nmax(Nmom,Ndat) ! bin number of last bin (computed from s_0)
      integer  :: i, j, k, l, kl, mn
      real(dp) :: S(Ndat), V1(Ndat), V1ERR(Ndat), A1(Ndat), A1ERR(Ndat)
      real(dp) :: COVVV(Ndat,Ndat), COVAA(Ndat,Ndat), COVAV(Ndat,Ndat)
      real(dp) :: sbin(Ndat), sfm2(2*Ndat), derr(2*Ndat)
      real(dp) :: errmat(2*Ndat+3,2*Ndat+3), jacobi(2*Ndat+3,2*Nsum)
      real(dp) :: mom_pi(Nsum), pifac, dpifac

!     Load rescaled OPAL data
      open(unit=9,&
          &file='data/udvec_specfunc_opal_hfagrescale_aug11.dat',&
          &status='old')
      do i = 1, Ndat
         read(9,*) S(i), V1(i), V1ERR(i)
      end do
      close(9)

      open(unit=9,&
          &file='data/udax_specfunc_opal_hfagrescale_aug11.dat',&
          &status='old')
      do i = 1, Ndat
         read(9,*) S(i), A1(i), A1ERR(i)
      end do
      close(9)

      open(unit=9,&
          &file='data/ud_covvv_covaa_covav_opal_hfagrescale_aug11.dat',&
          &status='old')
      do i = 1, Ndat
         do j = 1, Ndat
            read(9,*) k, l, COVVV(k,l), COVAA(k,l), COVAV(k,l)
         end do
      end do
      close(9)

!     Transform rescaled Opal data to Aleph normalisation
!     Additional factor 2Pi^2 in Kim's normalisation compared to Opal
      do i = 1, Ndat
         sbin(i) = S(i)
         sfm2(i)      = 12.d0*pi**2*SEW*Vud**2*Be*3.2d-2/stau*&
                       &real(wr00(dcmplx(sbin(i)/stau,0)),dp)*V1(i)
         sfm2(i+Ndat) = 12.d0*pi**2*SEW*Vud**2*Be*3.2d-2/stau*&
                       &real(wr00(dcmplx(sbin(i)/stau,0)),dp)*A1(i)
      end do

!     Renormalise error
      do i = 1, Ndat
         if (V1(i) /= 0.d0) then
            derr(i) = sfm2(i)/V1(i)*V1ERR(i)
         end if
         if (A1(i) /= 0.d0) then
            derr(i+Ndat) = sfm2(i+Ndat)/A1(i)*A1ERR(i)
         end if
      end do

!     Compute bin number corresponding to last bin
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            Nmax(k,l) = ceiling(31.25d0*s0s(k,l))
         end do
      end do

!     Calculate vector moments integrated up to s_0
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               mom(kl) = mom(kl) + stau/s0s(k,l)/Be*sfm2(i)*&
                        &real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))/&
                             &wr00(dcmplx(sbin(i)/stau,0)))
            end do
         end do
      end do

!     Get pion-pole contribution
      call pimom(Nsum,Ns0s,s0s,mom_pi,pifac,dpifac)

!     Calculate axialvector moments integrated up to s_0
!     Now kl starts with Nsum!
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               mom(kl) = mom(kl) + stau/s0s(k,l)/Be*sfm2(i+Ndat)*&
                        &real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))/&
                             &wr00(dcmplx(sbin(i)/stau,0)))
            end do
            mom(kl) = mom(kl) + mom_pi(kl-Nsum)
         end do
      end do

!     Initialise error matrix
      do i = 1, 2*Ndat+3
         do j = 1, 2*Ndat+3
            errmat(i,j) = 0.d0
         end do
      end do

!     Construct error matrix
      do i = 1, Ndat
         do j = 1, Ndat
            if ((V1(i)/=0.d0).and.(V1(j)/=0.d0)) then
               errmat(i,j) = sfm2(i)/V1(i)*sfm2(j)/V1(j)*COVVV(i,j)
            end if
!     Comment out to remove V/A cross correlations
            if ((A1(i)/=0.d0).and.(V1(j)/=0.d0)) then
               errmat(i+Ndat+1,j) = sfm2(i+Ndat)/A1(i)*sfm2(j)/V1(j)*&
                                   &COVAV(i,j)
            end if
            if ((A1(j)/=0.d0).and.(V1(i)/=0.d0)) then
               errmat(i,j+Ndat+1) = sfm2(j+Ndat)/A1(j)*sfm2(i)/V1(i)*&
                                   &COVAV(j,i)
            end if
!     Up to here
            if ((A1(i)/=0.d0).and.(A1(j)/=0.d0)) then
               errmat(i+Ndat+1,j+Ndat+1) = sfm2(i+Ndat)/A1(i)*&
                                          &sfm2(j+Ndat)/A1(j)*COVAA(i,j)
            end if
!     Switch off experimental correlations
!           if (i /= j) then
!              errmat(i,j) = 0.d0
!              errmat(i+Ndat+1,j+Ndat+1) = 0.d0
!           end if
!           errmat(i+Ndat+1,j) = 0.d0
!           errmat(i,j+Ndat+1) = 0.d0
         end do
      end do
      errmat(  Ndat+1,  Ndat+1) = dBe**2
      errmat(2*Ndat+2,2*Ndat+2) = dBe**2
      errmat(2*Ndat+3,2*Ndat+3) = dpifac**2

!     Initialise Jacobi matrix
      do i = 1, 2*Ndat+3
         do j = 1, 2*Nsum
            jacobi(i,j) = 0.d0
         end do
      end do

!     Construct Jacobi matrix; vector part
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               jacobi(i,kl) = stau/(s0s(k,l)*Be)*&
                             &real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))&
                                 &/wr00(dcmplx(sbin(i)/stau,0)))
            end do
            jacobi(Ndat+1,kl) = - mom(kl)/Be
         end do
      end do

!     Construct Jacobi matrix; axialvector part
!     Now kl starts with Nsum!
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               jacobi(i+Ndat+1,kl) = stau/(s0s(k,l)*Be)*&
                              &real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))&
                                  &/wr00(dcmplx(sbin(i)/stau,0)))
            end do
            jacobi(2*Ndat+2,kl) = (mom_pi(kl-Nsum) - mom(kl))/Be
            jacobi(2*Ndat+3,kl) =  mom_pi(kl-Nsum)/pifac
         end do
      end do

!     Calculate covariance matrix of relative moments
      do kl = 1, 2*Nsum
         do mn = 1, 2*Nsum
            do i = 1, 2*Ndat+3
               do j = 1, 2*Ndat+3
                  covmom(kl,mn) = covmom(kl,mn) +&
                 &jacobi(i,kl)*errmat(i,j)*jacobi(j,mn) 
               end do
            end do
         end do
      end do

      return
      end subroutine momex_va
!-----------------------------------------------------------------------
