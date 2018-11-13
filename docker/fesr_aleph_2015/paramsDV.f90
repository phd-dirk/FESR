      ! Fit results for DV parameters
      ! to be added as additional constraints in the chi^2
      module params_DV
         use num_const

         implicit none
!     Fit result for vector DV parameters
         real(dp), dimension(4) :: dvVval, dvVerr
         data dvVval(:) /3.46d-2, 0.691d0, 0.715d0, 2.640d0/
         data dvVerr(:) /1.07d-2, 0.179d0, 0.520d0, 0.284d0/

!  6    kapV       0.34609E-01   0.10682E-01   0.21004E-06   0.38481E-05
!  7    gamV       0.69137       0.17938       0.18388E-04   0.10989E-05
!  8    alpV       0.71512       0.51998       0.38967E-04  -0.47296E-06
!  9    betV        2.6403       0.28352       0.21887E-04   0.53870E-06

         real(dp), dimension(4,4) :: dvVcor
         data dvVcor(:,:) / 1.000d0,  0.987d0,  0.435d0, -0.401d0,&
                          & 0.987d0,  1.000d0,  0.372d0, -0.341d0,&
                          & 0.435d0,  0.372d0,  1.000d0, -0.993d0,&
                          &-0.401d0, -0.341d0, -0.993d0,  1.000d0/

! PARAMETER  CORRELATION COEFFICIENTS
!      NO.  GLOBAL     6     7     8     9
!       6  0.99042  1.000 0.987 0.435-0.401
!       7  0.98958  0.987 1.000 0.372-0.341
!       8  0.99444  0.435 0.372 1.000-0.993
!       9  0.99410 -0.401-0.341-0.993 1.000


!     Fit result for axialvector DV parameters
         real(dp), dimension(4) :: dvAval, dvAerr
         data dvAval(:) /0.2004d0, 1.480d0, -2.306d0, -1.989d0/
         data dvAerr(:) /0.0806d0, 0.181d0,  0.636d0,  0.330d0/

!  6    kapA       0.20039       0.80581E-01   0.12979E-04   0.32163E-06
!  7    gamA        1.4804       0.18124       0.39919E-04   0.90390E-06
!  8    alpA       -2.3063       0.63618       0.30810E-04   0.55565E-06
!  9    betA       -1.9885       0.33013       0.17444E-04   0.13152E-05

         real(dp), dimension(4,4) :: dvAcor
         data dvAcor(:,:) / 1.000d0,  0.980d0, -0.831d0,  0.804d0,&
                          & 0.980d0,  1.000d0, -0.742d0,  0.713d0,&
                          &-0.831d0, -0.742d0,  1.000d0, -0.998d0,&
                          & 0.804d0,  0.713d0, -0.998d0,  1.000d0/

! PARAMETER  CORRELATION COEFFICIENTS
!      NO.  GLOBAL     6     7     8     9
!       6  0.99466  1.000 0.980-0.831 0.804
!       7  0.99000  0.980 1.000-0.742 0.713
!       8  0.99952 -0.831-0.742 1.000-0.998
!       9  0.99943  0.804 0.713-0.998 1.000

      end module params_DV
