!     Definition of constants for the Renormalon Pole Model (RPM)
!     The following assumes number of flavours n_f = 3!
!     Dependence of residues on c_51 implemented

      module RPM_const
         use num_const

         implicit none
!     Residua for renormalon poles and polynomial coefficients
!        real(dp) :: cUV(3), cIR(5), cPO(3)
!     Main model with 1 UV and 2 IR poles and 2 polynomial terms
!        data cUV(1) /-0.01559637506686463d0/
!        data cIR(2:3) /3.157694392224263d0, -13.530201013824815d0/
!        data cPO(1:2) /7.808377283089274d-1, 7.658657314624812d-3/

!     Anomalous dimensions of D=2p operators (Only implemented for <aGG>)
         real(dp) :: cc_1(-3:6), cc_2(-3:6)
         data cc_1(2), cc_2(2) /-0.6111111111111111d0, 0.d0/

      contains

!-----------------------------------------------------------------------

!     Residues of the renormalon poles for the central RPM
!     as a function of c_51

         function cUV(n,c51)
            use num_const
            implicit none
            integer, intent(in)  :: n
            real(dp), intent(in) :: c51
            real(dp)             :: cUV

            cUV = 0.d0

            select case (n)
               case (1)
               cUV = -6.41040536064710827065d-2 +&
                     &1.71405224521570751843d-4*c51
            end select

         end function cUV

         function cIR(n,c51)
            use num_const
            implicit none
            integer, intent(in)  :: n
            real(dp), intent(in) :: c51
            real(dp)             :: cIR

            cIR = 0.d0

            select case (n)
               case (2)
               cIR =  2.71924081074743794779d0 +&
                     &1.54930594161429502669d-3*c51
               case (3)
               cIR = -3.52595444085400301258d1 +&
                     &7.67821321368025774151d-2*c51
            end select

         end function cIR

         function cPO(n,c51)
            use num_const
            implicit none
            integer, intent(in)  :: n
            real(dp), intent(in) :: c51
            real(dp)             :: cPO

            cPO = 0.d0

            select case (n)
               case (1)
               cPO = -1.22780380655402919290d0 +&
                     &7.09767326806695480791d-3*c51
               case (2)
               cPO = -0.74613717123668823940d0 +&
                     &2.66358950018128885880d-3*c51
            end select

         end function cPO

!-----------------------------------------------------------------------

         function gaUV(k,p)
            use num_const
            use rge_const

            implicit none
            integer, intent(in) :: k, p
            real(dp)            :: gaUV

            gaUV = k - 2.d0*p*be_2(3)/be_1(3)**2

         end function gaUV

!-----------------------------------------------------------------------

         function gaIR(k,p)
            use num_const
            use rge_const

            implicit none
            integer, intent(in) :: k, p
            real(dp)            :: gaIR

            gaIR = k + 2.d0*p*be_2(3)/be_1(3)**2

         end function gaIR

!-----------------------------------------------------------------------

         function bb_1(p)
            use num_const
            use rge_const

            implicit none
            integer, intent(in) :: p
            real(dp)            :: bb_1

            bb_1 = 2.d0*p/be_1(3)**3*(be_2(3)**2-be_1(3)*be_3(3))

         end function bb_1

!-----------------------------------------------------------------------

         function bb_2(p)
            use num_const
            use rge_const

            implicit none
            integer, intent(in) :: p
            real(dp)            :: bb_2

            real(dp) :: bb_1

            bb_1 = 2.d0*p/be_1(3)**3*(be_2(3)**2-be_1(3)*be_3(3))
            bb_2 = bb_1**2/2.d0-p/be_1(3)**4*(be_2(3)**3-2.d0*&
                  &be_1(3)*be_2(3)*be_3(3)+be_1(3)**2*be_4(3))

         end function bb_2
            
      end module RPM_const
