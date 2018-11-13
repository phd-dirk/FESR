
!     Compute experimental moments with correlations from OPAL data
!     Last change: Matthias Jamin, 26.3.2011

!-----------------------------------------------------------------------

!     Moment integrals of V/A spectral function

      subroutine momex_va(Nsum,Ns0s,s0s,mom,covmom)
      use params
      use opal98
      use weights_rho
      use weights_used

      implicit none
      integer,  intent(in) :: Nsum ! Total number of s_0's
      integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
      real(dp), dimension(Nmom,0:140) :: s0s ! Sets of used s_0's
      real(dp), intent(out) :: mom(2*Nsum), covmom(2*Nsum,2*Nsum)

      integer, parameter :: Ndat = 100 ! Number of Opal data points
      integer  :: Nmax(Nmom,Ndat) ! bin number of last bin (computed from s_0)
      integer  :: i, j, k, l, kl, mn
      real(dp) :: stauOpal = 1.777d0**2 ! M_tau^2(Opal)
      real(dp) :: sbin(Ndat), sfm2(2*Ndat), derr(2*Ndat),&
                 &corerr(2*Ndat,2*Ndat)
      real(dp) :: errmat(2*Ndat+3,2*Ndat+3), jacobi(2*Ndat+3,2*Nsum)
      real(dp) :: mom_pi(Nsum), pifac, dpifac

!     Remove error component of |V_ud|^2 and B_e
      do i = 1, Ndat
!        Relative pure data error
         derr(i)      = sqrt(V1ERR(i)**2-V1(i)**2*((0.08d0/17.83d0)**2+&
                            &(0.0008d0/0.9512d0)**2))
         derr(i+Ndat) = sqrt(A1ERR(i)**2-A1(i)**2*((0.08d0/17.83d0)**2+&
                            &(0.0008d0/0.9512d0)**2))
      end do
      
!     Transform Opal data to Aleph normalisation
      do i = 1, Ndat
         sbin(i) = S(i)
         sfm2(i)      = 6.d0*1.0194d0*0.9512d0*17.83d0*3.2d-2/stauOpal*&
                       &wr00(dcmplx(sbin(i)/stauOpal,0))*V1(i)
         sfm2(i+Ndat) = 6.d0*1.0194d0*0.9512d0*17.83d0*3.2d-2/stauOpal*&
                       &wr00(dcmplx(sbin(i)/stauOpal,0))*A1(i)
      end do

!     Renormalise error
      do i = 1, Ndat
         if (V1(i) /= 0.d0) then
            derr(i) = sfm2(i)/V1(i)*derr(i)
         end if
         if (A1(i) /= 0.d0) then
            derr(i+Ndat) = sfm2(i+Ndat)/A1(i)*derr(i+Ndat)
         end if
      end do

!     Construct correlation matrix
      do i = 1, Ndat
         do j = 1, Ndat
            corerr(i     ,j     ) = (CORRVV(i,j)+CORRVV(j,i))/2.d0
            corerr(i     ,j+Ndat) =  CORRVA(j,i)
!           corerr(i     ,j+Ndat) =  0.d0  ! Switch off cross-correlation
            corerr(i+Ndat,j     ) =  CORRVA(i,j)
!           corerr(i+Ndat,j     ) =  0.d0  ! Switch off cross-correlation
            corerr(i+Ndat,j+Ndat) = (CORRAA(i,j)+CORRAA(j,i))/2.d0
         end do
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
            errmat(i,j) = corerr(i,j)*derr(i)*derr(j)/1.d2
            errmat(i+Ndat+1,j) = corerr(i+Ndat,j)*derr(i+Ndat)*&
                                                 &derr(j)/1.d2
            errmat(i,j+Ndat+1) = corerr(i,j+Ndat)*derr(i)*&
                                                 &derr(j+Ndat)/1.d2
            errmat(i+Ndat+1,j+Ndat+1) = corerr(i+Ndat,j+Ndat)*&
                                       &derr(i+Ndat)*derr(j+Ndat)/1.d2
!     Switch off correlations
!           if (i /= j) then
!              errmat(i,j) = 0.d0
!           end if
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
