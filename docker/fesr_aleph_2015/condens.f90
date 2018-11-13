!     Definition of numerical input parameters
      module condens
         use params
         use rge_const

         implicit none
         real(dp), parameter :: atauBJ = astauBJ/pi ! Should be generalised!

!     Quark and gluon condensates at M_tau
         real(dp), parameter :: uumtau = -(0.272d0)**3 ! <uu>(M_tau)
         real(dp), parameter :: ddmtau = -(0.272d0)**3 ! <dd>(M_tau)
         real(dp), parameter :: kappa = 0.8d0          ! <ss>/<qq>
         real(dp), parameter :: ssmtau = kappa*uumtau  ! <ss>(M_tau)
         real(dp), parameter :: aGGmtau = 0.012d0      ! <aGG>(M_tau)

         real(dp), dimension(3), parameter :: qqmtau = (/&
                                             &uumtau, ddmtau, ssmtau /)

!     Invariant quark and gluon condensates

         real(dp), dimension(3), parameter :: qqinv = (/&
        &uumtau+3.d0*mumtau**3/(7.d0*pi**2)*(1.d0/atauBJ-53.d0/24.d0),&
        &ddmtau+3.d0*mdmtau**3/(7.d0*pi**2)*(1.d0/atauBJ-53.d0/24.d0),&
        &ssmtau+3.d0*msmtau**3/(7.d0*pi**2)*(1.d0/atauBJ-53.d0/24.d0) /)

!     For <aGG>_inv prepare beta- and gamma-functions

         real(dp), parameter :: beta3mtau = be_1(3)*atauBJ +&
        &be_2(3)*atauBJ**2 + be_3(3)*atauBJ**3 + be_4(3)*atauBJ**4
         real(dp), parameter :: gamm3mtau = ga_1(3)*atauBJ +&
        &ga_2(3)*atauBJ**2 + ga_3(3)*atauBJ**3 + ga_4(3)*atauBJ**4
         real(dp), parameter :: gam03mtau = - 2.d0 - 8.d0/3.d0*atauBJ
         
         real(dp), parameter :: aGGinv = beta3mtau/be_1(3)/&
                  &atauBJ*aGGmtau - 4.d0*gamm3mtau/be_1(3)*&
                  &(mumtau + mdmtau + msmtau*kappa)*uumtau +&
                  &3.d0/(4.d0*pi**2)*gam03mtau/be_1(3)*&
                  &(mumtau**4+mdmtau**4+msmtau**4)

      end module condens
