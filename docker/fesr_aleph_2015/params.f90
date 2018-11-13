!     Definition of numerical input parameters
      module params
         use num_const

         implicit none
!     Some particle masses
         real(dp), parameter :: mpim = 0.13957018d0 ! M_pi^-
         real(dp), parameter :: mkm =  0.493677d0   ! M_K^-
         real(dp), parameter :: mtau = 1.77686d0,&  ! PDG 2018
!        real(dp), parameter :: mtau = 1.77677d0,&  ! HFAG 2010
!        real(dp), parameter :: mtau = 1.77684d0,&  ! PDG 2008
!        real(dp), parameter :: mtau = 1.7769d0,&   ! Davier et al. 2008
!        real(dp), parameter :: mtau = 1.77703d0,&  ! ALEPH report 2005
!        real(dp), parameter :: mtau = 1.7770d0,&   ! OPAL 1998
!        real(dp), parameter :: mtau = 1.7765d0,&   ! to reproduce Tab. 23
                               &stau = mtau**2

!      Pseudoscalar resonance parameters
         real(dp), parameter ::  fpi = 92.21d-3,& ! PDG 2010
!        real(dp), parameter ::  fpi = 92.056d-3,&! Davier et al. 2008
!        real(dp), parameter ::  fpi = 92.10d-3,& ! to reproduce Tab. 23
!        real(dp), parameter ::  fpi = 94.00d-3,& ! OPAL 1998
                               &dfpi =  0.14d-3   ! PDG 2010
!                              &dfpi =  0.44d-3   ! Davier et al. 2008
!                              &dfpi =  0.44d-3   ! to reproduce Tab. 23
         real(dp), parameter ::  fk  = 1.198d0*fpi ! Kaon decay constant
!      Exited resonance parameters
         real(dp) ::  f1p = 2.20d-3, m1p = 1.300d0, g1p = 0.400d0
         real(dp) ::  f2p = 0.19d-3, m2p = 1.800d0, g2p = 0.210d0
         real(dp) ::  f1k = 21.4d-3, m1k = 1.460d0, g1k = 0.260d0
         real(dp) ::  f2k = 4.50d-3, m2k = 1.830d0, g2k = 0.250d0

!     Other parameters
         real(dp), parameter ::  Be = 17.815d0,& ! HFAG 2017
                               &dBe =  0.023d0
!        real(dp), parameter ::  Be = 17.827d0,& ! HFAG 2011
!                              &dBe =  0.000d0
!                              &dBe =  0.040d0
!        real(dp), parameter ::  Be = 17.85d0,&  ! PDG 2010
!                              &dBe =  0.05d0
!        real(dp), parameter ::  Be = 17.833d0,& ! Banerjee 2008
!                              &dBe =  0.030d0
!        real(dp), parameter ::  Be = 17.818d0,& ! Davier et al. 2008
!                              &dBe =  0.032d0
!        real(dp), parameter ::  Be = 17.810d0,& ! ALEPH report 2005
!                              &dBe =  0.039d0
!        real(dp), parameter ::  Be = 17.830d0,& ! OPAL 1998
!                              &dBe =  0.080d0
         real(dp), parameter ::&
!                               &RtauVex = 1.783d0,& ! Davier et al. 2008
!                              &dRtauVex = 0.011d0,&
!                               &RtauAex = 1.695d0,&
!                              &dRtauAex = 0.011d0,&
!                               &RtauVpAex = 3.479d0,&
!                              &dRtauVpAex = 0.011d0,&
!                               &RtauVpAex = 3.482d0,& ! ALEPH report 2005
!                              &dRtauVpAex = 0.014d0,&
!                               &RtauVpAex = 3.484d0,& ! OPAL 1998
!                              &dRtauVpAex = 0.024d0,&

!                               &RtauVpAex = 1.0d0,&
!                              &dRtauVpAex = 0.0d0,&
                                &RtauVmAex = 1.0d0,&
                               &dRtauVmAex = 0.0d0,&
                                &RtauVex = 1.0d0,&
                               &dRtauVex = 0.0d0,&
                                &RtauAex = 1.0d0,&
                               &dRtauAex = 0.0d0

         real(dp), parameter ::  Vud = 0.97420d0,& ! Towner, Hardy 2018
                               &dVud = 0.00021d0
!        real(dp), parameter ::  Vud = 0.97425d0,& ! Towner, Hardy 2009
!                              &dVud = 0.00022d0
!        real(dp), parameter ::  Vud = 0.97418d0,& ! CKMfitter 2008
!                              &dVud = 0.00019d0
!        real(dp), parameter ::  Vud = 0.9746d0,&  ! ALEPH report 2005
!                              &dVud = 0.0006d0
!        real(dp), parameter ::  Vud = 0.9753d0,&  ! OPAL 1998
!                              &dVud = 0.0006d0
!        real(dp), parameter ::  Vud = 0.9752d0,&  ! Andreas thesis
!                              &dVud = 0.0007d0
         real(dp), parameter ::  Vus = sqrt(1.d0-Vud**2) ! Unitarity

         real(dp), parameter ::  SEW = 1.0198d0,& ! EW radiative corr.
!        real(dp), parameter ::  SEW = 1.0194d0,& ! OPAL 1998
!        real(dp), parameter ::  SEW = 1.0201d0,& ! Value in V_us analysis
                               &dSEW = 0.0006d0,&
!                              &dSEW = 0.0040d0,& ! Andreas thesis
                            &deltaEW = 0.0010d0

!     Different input values for alpha_s(M_tau)
         real(dp), parameter :: astauMZ1180 = 0.3153862258436091d0,&
                               &astauMZ1184 = 0.3186045439284738d0,&
                               &astauMZ1190 = 0.3235272367246068d0,&
                               &astauBJ = 0.3156d0

!     Lambda_MSbar corresponding to alpha_s(M_Z) = 0.118
         real(dp), dimension(3:6), parameter ::&
            &lambda = (/ 0.3317858916330567d0,&
                        &0.2888660556427215d0,&
                        &0.2083638452037246d0,&
                        &0.2083638452037246d0 /) ! CORRECT VALUE MISSING!

!     Light quark masses
         real(dp), parameter :: mumtau = 2.8d-3,& ! m_up(M_tau)
                               &mdmtau = 5.0d-3,& ! m_down(M_tau)
                               &msmtau = 97.d-3   ! m_strange(M_tau)
!                              &msmtau = 98.4d-3  ! m_strange(M_tau)
!     Array notation for later use of quark masses
         real(dp), dimension(3), parameter ::&
                               &mq = (/ mumtau, mdmtau, msmtau/)

!     Uncertainties in c_41, c_51 and renormalisation scale
         real(dp), parameter :: dc41 =  25.d0,&   ! Aleph report 2005
!        real(dp), parameter :: dc41 =  50.d0,&   ! Andreas thesis
                               &dc51 = 283.d0,&
!                              &dc51 = 378.d0,&   ! Davier et al. 2008
                               &dmu2 =   1.0d0
!                              &dmu2 =   2.0d0    ! Davier et al. 2008
!                              &dmu2 =   2.0d0    ! Aleph report 2005
!                              &dmu2 =   1.5d0**2 ! Andreas thesis ?

      end module params
