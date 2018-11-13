
!     Compute experimental moments with correlations from ALEPH data
!     Detailed documentation in VAmomAleph.pdf
!     Last change: Matthias Jamin, 19.6.2014

!-----------------------------------------------------------------------

!     Compute moment integrals of V spectral function with
!     respect to an external value of R_tau,V

      subroutine momex_v(mom,covmom)
      use params
      use aleph_vv
      use weights_rho
      use weights_D
      use weights_used
      use s0s_fit

      implicit none
      real(dp) :: mom(Nsum), covmom(Nsum,Nsum)

      integer, parameter :: Ndat = 80 ! Number of Aleph 2014 data points
      integer  :: Nmax(Nmom,Ndat) ! maximal number of bins (computed from s_0)
      integer  :: i, j, k, l, kl, mn
      real(dp) :: sbin(Ndat), dsbin(Ndat), sfm2(Ndat), derr(Ndat)
      real(dp) :: corerr(Ndat,Ndat), errmat(Ndat+2,Ndat+2)
      real(dp) :: jacobi(Ndat+2,Nsum), wratio(Ndat,Nsum)

!     Load ALEPH 2014 vector data
      call aleph_v(sbin,dsbin,sfm2,derr,corerr)

!     Initialise used set of s_0's
      call init_s0s_fit(Nsum,Ns0s,s0s)

!     Set Nmax(k,l) to the bin number of the last included bin
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            i = Ndat
            do
               i = i - 1
               if (abs(s0s(k,l)-1.d-6-sbin(i)) < dsbin(i)/2.d0) exit
            end do
               Nmax(k,l) = i
!              write (*,*) Nmax(k,l)
         end do
      end do

!     Calculate relative moments integrated up to s_0
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
!     Ratio of moments at centre of bins
!              wratio(i,kl) = real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))/&
!                            &wr00(dcmplx(sbin(i)/stau,0)))
!     Ratio of integrated moments
               wratio(i,kl) = s0s(k,l)/stau*real(&
               &(wgtD(k,dcmplx((sbin(i)-dsbin(i)/2.d0)/s0s(k,l),0)) -&
                &wgtD(k,dcmplx((sbin(i)+dsbin(i)/2.d0)/s0s(k,l),0)))/&
               &(wD00(dcmplx((sbin(i)-dsbin(i)/2.d0)/stau,0)) -&
                &wD00(dcmplx((sbin(i)+dsbin(i)/2.d0)/stau,0))))

               mom(kl) = mom(kl) + stau/s0s(k,l)/Be*sfm2(i)*wratio(i,kl)
            end do
            if (k /= 1) then
               mom(kl) = mom(kl)/RtauVex
            end if
         end do
      end do

!     Initialise error matrix
      do i = 1, Ndat+2
         do j = 1, Ndat+2
            errmat(i,j) = 0.d0
         end do
      end do

!     Construct error matrix
      do i = 1, Ndat
         do j = 1, Ndat
            errmat(i,j) = corerr(i,j)*derr(i)*derr(j)/1.d2
!     Switch off correlations
!           if (i /= j) then
!              errmat(i,j) = 0.d0
!           end if
         end do
      end do
      errmat(Ndat+1,Ndat+1) = dBe**2
      errmat(Ndat+2,Ndat+2) = dRtauVex**2

!     Initialise Jacobi matrix
      do i = 1, Ndat+2
         do j = 1, Nsum
            jacobi(i,j) = 0.d0
         end do
      end do

!     Construct Jacobi matrix
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               jacobi(i,kl) = stau/(s0s(k,l)*Be)*wratio(i,kl)
               if (k /= 1) then
                  jacobi(i,kl) = jacobi(i,kl)/RtauVex
               end if
            end do
            jacobi(Ndat+1,kl) = - mom(kl)/Be
            if (k /= 1) then
               jacobi(Ndat+2,kl) = - mom(kl)/RtauVex
            end if
         end do
      end do

!     Write Jacobi matrix
!     do i = 1, Ndat+2
!        do j = 1, Nsum
!           write(*,*) jacobi(i,j)
!        end do
!     end do

!     Calculate covariance matrix of relative moments
      do kl = 1, Nsum
         do mn = 1, Nsum
            do i = 1, Ndat+2
               do j = 1, Ndat+2
                  covmom(kl,mn) = covmom(kl,mn) +&
                 &jacobi(i,kl)*errmat(i,j)*jacobi(j,mn) 
               end do
            end do
         end do
      end do

      return
      end subroutine momex_v
!-----------------------------------------------------------------------

!     Compute moment integrals of V+A spectral function

      subroutine momex_vpa(mom,covmom)
      use params
      use aleph_vpa
      use weights_rho
      use weights_D
      use weights_used
      use s0s_fit

      implicit none
      real(dp) :: mom(Nsum), covmom(Nsum,Nsum)

      integer, parameter :: Ndat = 80 ! Number of Aleph 2014 data points
      integer  :: Nmax(Nmom,Ndat) ! maximal number of bins (computed from s_0)
      integer  :: i, j, k, l, kl, mn
      real(dp) :: sbin(Ndat), dsbin(Ndat), sfm2(Ndat), derr(Ndat)
      real(dp) :: corerr(Ndat,Ndat), errmat(Ndat+2,Ndat+2)
      real(dp) :: jacobi(Ndat+2,Nsum), wratio(Ndat,Nsum)
      real(dp) :: mom_pi(Nsum), pifac, dpifac

!     Load ALEPH 2014 vector plus axialvector data
      call aleph_vplusa(sbin,dsbin,sfm2,derr,corerr)

!     Some renormalisation of R_tau,V+A
      do i = 1, Ndat
         sfm2(i) = 0.99743669d0*sfm2(i)
         derr(i) = 0.99743669d0*derr(i)
      end do

!     Initialise used set of s_0's
      call init_s0s_fit(Nsum,Ns0s,s0s)

!     Set Nmax(k,l) to the bin number of the last included bin
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            i = Ndat
            do
               i = i - 1
               if (abs(s0s(k,l)-1.d-6-sbin(i)) < dsbin(i)/2.d0) exit
            end do
               Nmax(k,l) = i
!              write (*,*) Nmax(k,l)
         end do
      end do

!     Get pion-pole contribution
      call pimom(Nsum,Ns0s,s0s,mom_pi,pifac,dpifac)
!     write (*,*) "pion_pole: ", mom_pi(1)

!     Calculate moments integrated up to s_0
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
!     Ratio of moments at centre of bins
!              wratio(i,kl) = real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))/&
!                            &wr00(dcmplx(sbin(i)/stau,0)))
!     Ratio of integrated moments
               wratio(i,kl) = s0s(k,l)/stau*real(&
               &(wgtD(k,dcmplx((sbin(i)-dsbin(i)/2.d0)/s0s(k,l),0)) -&
                &wgtD(k,dcmplx((sbin(i)+dsbin(i)/2.d0)/s0s(k,l),0)))/&
               &(wD00(dcmplx((sbin(i)-dsbin(i)/2.d0)/stau,0)) -&
                &wD00(dcmplx((sbin(i)+dsbin(i)/2.d0)/stau,0))))

               mom(kl) = mom(kl) + stau/s0s(k,l)/Be*sfm2(i)*wratio(i,kl)
            end do
            mom(kl) = mom(kl) + mom_pi(kl)
         end do
      end do

!     Initialise error matrix
      do i = 1, Ndat+2
         do j = 1, Ndat+2
            errmat(i,j) = 0.d0
         end do
      end do

!     Construct error matrix
      do i = 1, Ndat
         do j = 1, Ndat
            errmat(i,j) = corerr(i,j)*derr(i)*derr(j)/1.d2
!     Switch off correlations
!           if (i /= j) then
!              errmat(i,j) = 0.d0
!           end if
         end do
      end do
      errmat(Ndat+1,Ndat+1) = dBe**2
      errmat(Ndat+2,Ndat+2) = dpifac**2

!     Initialise Jacobi matrix
      do i = 1, Ndat+2
         do j = 1, Nsum
            jacobi(i,j) = 0.d0
         end do
      end do


!     Construct Jacobi matrix
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               jacobi(i,kl) = stau/(s0s(k,l)*Be)*wratio(i,kl)
            end do
            jacobi(Ndat+1,kl) = (mom_pi(kl) - mom(kl))/Be
            jacobi(Ndat+2,kl) =  mom_pi(kl)/pifac
         end do
      end do

!     Calculate covariance matrix of relative moments
      do kl = 1, Nsum
         do mn = 1, Nsum
            do i = 1, Ndat+2
               do j = 1, Ndat+2
                  covmom(kl,mn) = covmom(kl,mn) +&
                 &jacobi(i,kl)*errmat(i,j)*jacobi(j,mn)
               end do
            end do
         end do
      end do

      return
      end subroutine momex_vpa
!-----------------------------------------------------------------------

!     Compute moment integrals of V-A spectral function

      subroutine momex_vma(mom,covmom)
      use params
      use aleph_vma
      use weights_rho
      use weights_D
      use weights_used
      use s0s_fit

      implicit none
      real(dp) :: mom(Nsum), covmom(Nsum,Nsum)

      integer, parameter :: Ndat = 80 ! Number of Aleph 2014 data points
      integer  :: Nmax(Nmom,Ndat) ! maximal number of bins (computed from s_0)
      integer  :: i, j, k, l, kl, mn
      real(dp) :: sbin(Ndat), dsbin(Ndat), sfm2(Ndat), derr(Ndat)
      real(dp) :: corerr(Ndat,Ndat), errmat(Ndat+2,Ndat+2)
      real(dp) :: jacobi(Ndat+2,Nsum), wratio(Ndat,Nsum)
      real(dp) :: mom_pi(Nsum), pifac, dpifac

!     Load ALEPH 2014 vector minus axialvector data
      call aleph_vminusa(sbin,dsbin,sfm2,derr,corerr)

!     Some renormalisation of R_tau,V-A
      do i = 1, Ndat
         sfm2(i) = 0.99363d0*sfm2(i)
         derr(i) = 0.99363d0*derr(i)
      end do

!     Initialise used set of s_0's
      call init_s0s_fit(Nsum,Ns0s,s0s)

!     Set Nmax(k,l) to the bin number of the last included bin
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            i = Ndat
            do
               i = i - 1
               if (abs(s0s(k,l)-1.d-6-sbin(i)) < dsbin(i)/2.d0) exit
            end do
               Nmax(k,l) = i
!              write (*,*) Nmax(k,l)
         end do
      end do

!     Get pion-pole contribution
      call pimom(Nsum,Ns0s,s0s,mom_pi,pifac,dpifac)
!     write (*,*) "pion_pole: ", mom_pi(1)

!     Calculate moments integrated up to s_0
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
!     Ratio of moments at centre of bins
!              wratio(i,kl) = real(wgtr(k,dcmplx(sbin(i)/s0s(k,l),0))/&
!                            &wr00(dcmplx(sbin(i)/stau,0)))
!     Ratio of integrated moments
               wratio(i,kl) = s0s(k,l)/stau*real(&
               &(wgtD(k,dcmplx((sbin(i)-dsbin(i)/2.d0)/s0s(k,l),0)) -&
                &wgtD(k,dcmplx((sbin(i)+dsbin(i)/2.d0)/s0s(k,l),0)))/&
               &(wD00(dcmplx((sbin(i)-dsbin(i)/2.d0)/stau,0)) -&
                &wD00(dcmplx((sbin(i)+dsbin(i)/2.d0)/stau,0))))

               mom(kl) = mom(kl) + stau/s0s(k,l)/Be*sfm2(i)*wratio(i,kl)
            end do
            mom(kl) = mom(kl) - mom_pi(kl)
         end do
      end do

!     Initialise error matrix
      do i = 1, Ndat+2
         do j = 1, Ndat+2
            errmat(i,j) = 0.d0
         end do
      end do

!     Construct error matrix
      do i = 1, Ndat
         do j = 1, Ndat
            errmat(i,j) = corerr(i,j)*derr(i)*derr(j)/1.d2
!     Switch off correlations
!           if (i /= j) then
!              errmat(i,j) = 0.d0
!           end if
         end do
      end do
      errmat(Ndat+1,Ndat+1) = dBe**2
      errmat(Ndat+2,Ndat+2) = dpifac**2

!     Initialise Jacobi matrix
      do i = 1, Ndat+2
         do j = 1, Nsum
            jacobi(i,j) = 0.d0
         end do
      end do


!     Construct Jacobi matrix
      kl = 0
      do k = 1, Nmom
         do l = 1, Ns0s(k)
            kl = kl + 1
            do i = 1, Nmax(k,l)
               jacobi(i,kl) = stau/(s0s(k,l)*Be)*wratio(i,kl)
            end do
            jacobi(Ndat+1,kl) = - (mom_pi(kl) + mom(kl))/Be
            jacobi(Ndat+2,kl) = - mom_pi(kl)/pifac
         end do
      end do

!     Calculate covariance matrix of relative moments
      do kl = 1, Nsum
         do mn = 1, Nsum
            do i = 1, Ndat+2
               do j = 1, Ndat+2
                  covmom(kl,mn) = covmom(kl,mn) +&
                 &jacobi(i,kl)*errmat(i,j)*jacobi(j,mn)
               end do
            end do
         end do
      end do

      return
      end subroutine momex_vma
!-----------------------------------------------------------------------
