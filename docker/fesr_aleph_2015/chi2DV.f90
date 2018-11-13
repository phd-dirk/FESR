
!     Calculate a chi^2 for the DV contributions
!     To be added to the chi^2 of the fit as an additional constraint
!     Last change: Matthias Jamin, 28.4.2010

!-----------------------------------------------------------------------

!     Chi^2 for DV for vector spectral function

      function chi2DVV(kaV,gaV,alV,beV)
      use params_DV
      use matrix

      implicit none
      real(dp), intent(in)  :: kaV, gaV, alV, beV
      real(dp) :: chi2DVV

      integer :: i, j, ierr
      real(dp), dimension(4)   :: dvdif
      real(dp), dimension(4,4) :: dverr, dvinv
      real(dp) :: sun

!     Set up paramter differences
      dvdif(1) = kaV - dvVval(1)
      dvdif(2) = gaV - dvVval(2)
      dvdif(3) = alV - dvVval(3)
      dvdif(4) = beV - dvVval(4)
      
!     Set up error matrix
      do i = 1, 4
         do j = 1, 4
            dverr(i,j) = dvVcor(i,j)*dvVerr(i)*dvVerr(j)
         end do
      end do

!     Invert error matrix
      call matinv(dverr,dvinv,4,ierr)

!     Compute the chi^2
      sun = 0.d0
      do i = 1, 4
         do j = 1, 4
            sun = sun + dvdif(i)*dvinv(i,j)*dvdif(j)
         end do
      end do

      chi2DVV = sun

      return
      end function chi2DVV
!-----------------------------------------------------------------------

!     Chi^2 for DV for axialvector spectral function

      function chi2DVA(kaA,gaA,alA,beA)
      use params_DV
      use matrix

      implicit none
      real(dp), intent(in)  :: kaA, gaA, alA, beA
      real(dp) :: chi2DVA

      integer :: i, j, ierr
      real(dp), dimension(4)   :: dvdif
      real(dp), dimension(4,4) :: dverr, dvinv
      real(dp) :: sun

!     Set up paramter differences
      dvdif(1) = kaA - dvAval(1)
      dvdif(2) = gaA - dvAval(2)
      dvdif(3) = alA - dvAval(3)
      dvdif(4) = beA - dvAval(4)
      
!     Set up error matrix
      do i = 1, 4
         do j = 1, 4
            dverr(i,j) = dvAcor(i,j)*dvAerr(i)*dvAerr(j)
         end do
      end do

!     Invert error matrix
      call matinv(dverr,dvinv,4,ierr)

!     Compute the chi^2
      sun = 0.d0
      do i = 1, 4
         do j = 1, 4
            sun = sun + dvdif(i)*dvinv(i,j)*dvdif(j)
         end do
      end do

      chi2DVA = sun

      return
      end function chi2DVA
!-----------------------------------------------------------------------
