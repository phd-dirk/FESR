!-----------------------------------------------------------------------

!     Write out experimental moments of ALEPH spectral functions

      subroutine writeVAmomAleph()
      use params
      use s0s_out
      use aleph_vpa
      use aleph_vma

      implicit none

      integer  :: i, j
      integer, parameter :: Ndat = 80 ! Number of Aleph 2014 data points
      real(dp) :: mom(Nsum), covmom(Nsum,Nsum)
      real(dp) :: cormat(Nsum,Nsum) ! Correlation matrix
      real(dp) :: sbin(Ndat), dsbin(Ndat), sfm2(Ndat), derr(Ndat)
      real(dp) :: corerr(Ndat,Ndat)
!     real(dp) :: mom(2*Nsum), covmom(2*Nsum,2*Nsum)
!     real(dp) :: cormat(2*Nsum,2*Nsum) ! Correlation matrix
!     real(dp) :: mom_pi(Nmom,Ns0s), pifac, dpifac

!     Initialise used set of s_0's
      call init_s0s_out(Nsum,Ns0s,s0s)

!     Load ALEPH 2014 vector plus axialvector data
      call aleph_vplusa(sbin,dsbin,sfm2,derr,corerr)
!     Load ALEPH 2014 vector minus axialvector data
!     call aleph_vminusa(sbin,dsbin,sfm2,derr,corerr)

!     write(*,*) ''
!     write(*,*) 'w(x) = 1'
!     write(*,FMT='(A8,5F10.5)') "s_0's = ",&
!                  &s0s(1,1), s0s(1,2), s0s(1,3), s0s(1,4), s0s(1,5)
!     write(*,*) 'w(x) = x'
!     write(*,FMT='(A8,3F10.5)') "s_0's = ",&
!                            &s0s(2,1), s0s(2,2), s0s(2,3)
!     write(*,*) 'w(x) = x^2'
!     write(*,FMT='(A8,2F10.5)') "s_0's = ",&
!                            &s0s(3,1), s0s(3,2)


!     Initialise moment vector and covariance matrix
      do i = 1, Nsum
         mom(i) = 0.d0
         do j = 1, Nsum
            covmom(i,j) = 0.d0
         end do
      end do


      write(*,*) ''
!     Compute moments and covariances for V, A, V+A, V-A spectral functions

!     A SPECTRAL MOMENTS
!     write(*,*) 'V SPECTRAL FUNCTIONS'
!     call    momex_v(mom,covmom)
!     call relmomex_v(Nmom,Ns0s,s0s,mom,covmom)

!     A SPECTRAL MOMENTS
!     write(*,*) 'A SPECTRAL FUNCTIONS'
!     call    momex_a(Nmom,Ns0s,s0s,mom,covmom)
!     call relmomex_a(Nmom,Ns0s,s0s,mom,covmom)

!     V+A SPECTRAL MOMENTS
!     write(*,*) 'V+A SPECTRAL FUNCTIONS'
      call    momex_vpa(mom,covmom)
!     call relmomex_vplusa(Nmom,Ns0s,s0s,mom,covmom)

!     V-A SPECTRAL MOMENTS
!     write(*,*) 'V-A SPECTRAL FUNCTIONS'
!     call    momex_vma(mom,covmom)
!     call relmomex_vminusa(Nmom,Ns0s,s0s,mom,covmom)


      do i = 1, Nsum
         do j = 1, Nsum
            cormat(i,j) = covmom(i,j)/sqrt(covmom(i,i)*covmom(j,j))
         end do
      end do

!     Write out results

!     write(*,*) ''
!     write(*,FMT='(A33,F8.5)')&
!           &'Spectral function moments: s_0 = ', s0s(l)
!     write(*,*) '   (0,0)     (1,0)     (1,1)     (1,2)     (1,3)&
!     &      w_M'
!     write(*,FMT='(5F14.8)') (mom(i), i=1,Nsum)
!     write(*,FMT='(5F14.8)') (sqrt(covmom(i,i)), i=1,Nsum)

!     write(*,*) ''
!     do i = 1, Nsum
!        write(*,FMT='(25F14.8)') (covmom(i,j), j=1,Nsum)
!     end do

!     write(*,*) ''
!     do i = 1, Nsum
!        write(*,FMT='(36F14.8)') (cormat(i,j), j=1,Nsum)
!     end do

      do i = 1, Nsum
         write(*,FMT='(78F14.8)') sqrt(s0s(1,i)), mom(i),&
                                 &sqrt(covmom(i,i))
      end do


!     Pion pole contribution

!     call pimom(Nmom,Ns0s,s0s,mom_pi,pifac,dpifac)

!     write(*,*) ''
!     write(*,*) ''
!     write(*,*) 'Separate pion-pole contribution'
!     write(*,*) '   (0,0)     (1,0)     (1,1)     (1,2)     (1,3)&
!     &      w_M'

!     do l = 1, Ns0s
!     write(*,*) ''
!     write(*,FMT='(6F10.5)') mom_pi(1,l), mom_pi(2,l), mom_pi(3,l),&
!                            &mom_pi(4,l), mom_pi(5,l), mom_pi(6,l)
!     write(*,FMT='(6F10.5)')&
!            &dpifac/pifac*mom_pi(1,l), dpifac/pifac*mom_pi(2,l),&
!            &dpifac/pifac*mom_pi(3,l), dpifac/pifac*mom_pi(4,l),&
!            &dpifac/pifac*mom_pi(5,l), dpifac/pifac*mom_pi(6,l)
!     end do

      return
      end subroutine writeVAmomAleph
!-----------------------------------------------------------------------
