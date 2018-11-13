!-----------------------------------------------------------------------

!     Write out experimental moments of OPAL spectral functions

      subroutine writeVAmomOpal(scheme,astau,aGGinv,rhoD6V,rhoD6A,&
                &c8D8V,c8D8A,cV0,kaV,gaV,alV,beV,kaA,gaA,alA,beA)
      use params
      use s0s_out

      implicit none

      integer  :: k, l, kl
      character(len=2)  :: scheme
      character(len=1)  :: k_char
      character(len=18) :: OutFileName
      real(dp) :: astau, aGGinv, rhoD6V, rhoD6A, c8D8V, c8D8A,&
                 &kaV, gaV, alV, beV, kaA, gaA, alA, beA
      real(dp), dimension(0:12,1:12,3:6) :: cV0
      real(dp) :: momex(2*Nsum), covex(2*Nsum,2*Nsum)
      real(dp) :: momthV(Nsum), momthA(Nsum)

!     Initialise used set of s_0's
      call init_s0s_out(Nsum,Ns0s,s0s)

!     Initialise moment vector and covariance matrix
      do k = 1, 2*Nsum
         momex(k) = 0.d0
         do l = 1, 2*Nsum
            covex(k,l) = 0.d0
         end do
      end do

!     OPAL MOMENTS
!     Compute experimental moments and covariances for V/A spectral function
      call momex_va(Nsum,Ns0s,s0s,momex,covex)

      kl = 0
!     Write experimental vector moments
      do k = 1, Nmom

         k_char = char(k+48)
         OutFileName = 'figs/momex_v_'// k_char //'.dat'
         Open (unit=9, file=OutFileName)

         do l = 1, Ns0s(k)
            kl = kl + 1
            write(9,FMT='(F6.3,2F12.6)') s0s(k,l), momex(kl),&
                                        &sqrt(covex(kl,kl))
         end do
         close (9) 
      end do

!     Write experimental axialvector moments
      do k = 1, Nmom

         k_char = char(k+48)
         OutFileName = 'figs/momex_a_'// k_char //'.dat'
         Open (unit=9, file=OutFileName)

         do l = 1, Ns0s(k)
            kl = kl + 1
            write(9,FMT='(3F11.6)') s0s(k,l), momex(kl),&
                                   &sqrt(covex(kl,kl))
         end do
         close (9)
      end do

!     Compute theoretical moments and covariances for V spectral function
      call momth_va(Nsum,Ns0s,s0s,scheme, 1,cV0,1.d0,astau,aGGinv,&
                   &rhoD6V,c8D8V,kaV,gaV,alV,beV,momthV)

      kl = 0
!     Write theoretical vector moments
      do k = 1, Nmom

         k_char = char(k+48)
         OutFileName = 'figs/momth_v_'// k_char //'.dat'
         Open (unit=9, file=OutFileName)

         do l = 1, Ns0s(k)
            kl = kl + 1
            write(9,FMT='(F6.3,F12.6)') s0s(k,l), momthV(kl)
         end do
         close (9)
      end do


!     Compute theoretical moments and covariances for A spectral function
      call momth_va(Nsum,Ns0s,s0s,scheme,-1,cV0,1.d0,astau,aGGinv,&
                   &rhoD6A,c8D8A,kaA,gaA,alA,beA,momthA)

      kl = 0
!     Write theoretical axialvector moments
      do k = 1, Nmom

         k_char = char(k+48)
         OutFileName = 'figs/momth_a_'// k_char //'.dat'
         Open (unit=9, file=OutFileName)

         do l = 1, Ns0s(k)
            kl = kl + 1
            write(9,FMT='(F6.3,F12.6)') s0s(k,l), momthA(kl)
         end do
         close (9)
      end do

      return
      end subroutine writeVAmomOpal
!-----------------------------------------------------------------------
