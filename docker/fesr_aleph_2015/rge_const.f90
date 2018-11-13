!     Definition of renormalisation group constants
!     Last change: Matthias Jamin, 5.4.2018
      module rge_const
         use num_const

         implicit none
         integer :: f  !  Number of flavours
!     beta-function coefficients beta_1 to beta_5
         real(dp), dimension(3:6), parameter ::&
            &be_1 = (/ (5.5d0 - f/3.d0, f=3,6) /),&
            &be_2 = (/ (12.75d0 - 19.d0/12.d0*f, f=3,6) /),&
            &be_3 = (/ (2857.d0/64.d0 - 5033.d0/576.d0*f +&
                       &325.d0/1728.d0*f*f, f=3,6) /),&
            &be_4 = (/ (149753.d0/768.d0 + 891.d0/32.d0*z3 -&
                       &(1078361.d0/20736.d0 + 1627.d0/864.d0*z3)*f +&
                       &(50065.d0/20736.d0 + 809.d0/1296.d0*z3)*f*f +&
                       &1093.d0/93312.d0*f**3, f=3,6) /),&
            &be_5 = (/ (8157455.d0/8192.d0 + 621885.d0/1024.d0*z3 -&
                       &88209.d0/1024.d0*z4 - 144045.d0/256.d0*z5 -&
                       &(336460813.d0/995328.d0 + 1202791.d0/10368.d0*&
                        &z3 - 33935.d0/3072.d0*z4 - 1358995.d0/&
                        &13824.d0*z5)*f +&
                       &(25960913.d0/995328.d0 + 698531.d0/41472.d0*z3&
                        &- 5263.d0/2304.d0*z4 - 5965.d0/648.d0*z5)*f*f&
                       &- (630559.d0/2985984.d0 + 24361.d0/62208.d0*z3&
                        &- 809.d0/6912.d0*z4 - 115.d0/1152.d0*z5)*f**3&
                       &+ (1205.d0/1492992.d0 - 19.d0/5184.d0*z3)*&
                        &f**4, f=3,6) /)

!     Auxiliary constants for 4-loop running of alpha_s
         real(dp), dimension(3:6), parameter ::&
            &be_ratio_21 = (/ (be_2(f)/be_1(f), f=3,6) /),&
            &be_ratio_31 = (/ (be_3(f)/be_1(f), f=3,6) /),&
            &be_ratio_41 = (/ (be_4(f)/be_1(f), f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &be_ratio_14 = (/ (be_1(f)/be_4(f), f=3,6) /),&
            &be_ratio_24 = (/ (be_2(f)/be_4(f), f=3,6) /),&
            &be_ratio_34 = (/ (be_3(f)/be_4(f), f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &pp = (/ (be_ratio_24(f)-be_ratio_34(f)**2/3.d0, f=3,6) /),&
            &qq = (/ (be_ratio_14(f)-be_ratio_24(f)*be_ratio_34(f)/&
                     &3.d0+2.d0*be_ratio_34(f)**3/27.d0, f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &dis = (/ ((qq(f)/2.d0)**2+(pp(f)/3.d0)**3, f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &qr_1 = (/ ( exp(log(sqrt(dis(f))-qq(f)/2.d0)/3.d0),&
                      & f=3,6) /),&
            &qr_2 = (/ (-exp(log(sqrt(dis(f))+qq(f)/2.d0)/3.d0),&
                      & f=3,6) /)

         real(dp), dimension(3:6) ::&
            &qw_1 = (/ (qr_1(f)+qr_2(f)+2.d0*be_ratio_34(f)/3.d0,&
                      &f=3,6) /),&
            &qw_2 = (/ (qr_1(f)**2+qr_2(f)**2-qr_1(f)*qr_2(f)+(qr_1(f)&
                       &+qr_2(f))*be_ratio_34(f)/3.d0+be_ratio_34(f)**2&
                       &/9.d0, f=3,6) /),&
            &qw_3 = (/ (be_ratio_34(f)/3.d0-qr_1(f)-qr_2(f), f=3,6) /)

!     Auxiliary constants for 5-loop running of alpha_s

         real(dp), dimension(3:6), parameter ::&
            &be_ratio_15 = (/ (be_1(f)/be_5(f), f=3,6) /),&  ! e
            &be_ratio_25 = (/ (be_2(f)/be_5(f), f=3,6) /),&  ! d
            &be_ratio_35 = (/ (be_3(f)/be_5(f), f=3,6) /),&  ! c
            &be_ratio_45 = (/ (be_4(f)/be_5(f), f=3,6) /)    ! b

         real(dp), dimension(3:6), parameter ::&
            &p5loop = (/ (be_ratio_35(f) -&
                         &0.375d0*be_ratio_45(f)**2, f=3,6) /),&
            &q5loop = (/ (be_ratio_45(f)**3/8.d0 - be_ratio_45(f)*&
                         &be_ratio_35(f)/2.d0 +&
                         &be_ratio_25(f), f=3,6) /),&
            &De0 = (/ (be_ratio_35(f)**2 - 3.d0*be_ratio_45(f)*&
                      &be_ratio_25(f)+12.d0*be_ratio_15(f), f=3,6) /),&
            &De1 = (/ (2.d0*be_ratio_35(f)**3 - 9.d0*be_ratio_45(f)*&
                      &be_ratio_35(f)*be_ratio_25(f) + 27.d0*&
                      &be_ratio_45(f)**2*be_ratio_15(f) + 27.d0*&
                      &be_ratio_25(f)**2 - 72.d0*be_ratio_35(f)*&
                      &be_ratio_15(f), f=3,6) /)

         complex(dc), dimension(3:6), parameter ::&
            &Q5l = (/ (((De1(f) + zsqrt(cmplx(De1(f)**2-4.d0*De0(f)**3,&
                       &0.d0,dc)))/2.d0)**(1.d0/3.d0), f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &S5l = (/ (real(zsqrt(cmplx(-2.d0/3.d0*p5loop(f),0.d0,dc) +&
                      &(Q5l(f)+De0(f)/Q5l(f))/3.d0))/2.d0, f=3,6) /)

!     roots of the quartic equation: Re(w), |w|^2
         real(dp), dimension(3:6), parameter ::&
            &Rew1 = (/ (- be_ratio_45(f)/4.d0 - S5l(f), f=3,6) /),&
            &Absw12 = (/ (be_ratio_45(f)**2/16.d0 + p5loop(f)/2.d0 -&
                         &q5loop(f)/4.d0/S5l(f) + be_ratio_45(f)*&
                         &S5l(f)/2.d0 + 2.D0*S5l(f)**2, f=3,6) /),&
            &Rew2 = (/ (- be_ratio_45(f)/4.d0 + S5l(f), f=3,6) /),&
            &Absw22 = (/ (be_ratio_45(f)**2/16.d0 + p5loop(f)/2.d0 +&
                         &q5loop(f)/4.d0/S5l(f) - be_ratio_45(f)*&
                         &S5l(f)/2.d0 + 2.D0*S5l(f)**2, f=3,6) /)

!     parameters in the partial fraction decomposition
         real(dp), dimension(3:6), parameter ::&
            &al1 = (/ (1.d0/(Absw12(f)*Absw22(f)), f=3,6) /),&
            &al2 = (/ (2.d0*(Rew1(f)/Absw12(f) + Rew2(f)/Absw22(f))/&
                      &(Absw12(f)*Absw22(f)), f=3,6) /),&
            &al3 = (/ ((Absw12(f)**2 - Absw12(f)*(Absw22(f) + 4.d0*&
                      &Rew1(f)*(3.d0*Rew1(f) - 2.d0*Rew2(f))) + 4.d0*&
                      &Rew1(f)**2*(Absw22(f) + 4.d0*Rew1(f)*(Rew1(f) -&
                      &Rew2(f))))/(Absw12(f)**2*(Absw12(f)**2 +&
                      &Absw22(f)*(Absw22(f) + 4.d0*Rew1(f)*(Rew1(f) -&
                      &Rew2(f))) - 2.d0*Absw12(f)*(Absw22(f) + 2.d0*&
                      &(Rew1(f) - Rew2(f))*Rew2(f)))), f=3,6) /),&
            &al4 = (/ ((-2.d0*Rew1(f)*(Absw22(f) + 4.d0*Rew1(f)*&
                      &(Rew1(f) - Rew2(f))) + Absw12(f)*(4.d0*Rew1(f) -&
                      &2.d0*Rew2(f)))/(Absw12(f)**2*(Absw12(f)**2 +&
                      &Absw22(f)*(Absw22(f) + 4.d0*Rew1(f)*(Rew1(f) -&
                      &Rew2(f))) - 2.d0*Absw12(f)*(Absw22(f) + 2.d0*&
                      &(Rew1(f) - Rew2(f))*Rew2(f)))), f=3,6) /),&
            &al5 = (/ ((Absw22(f)*(-Absw12(f) + Absw22(f)) + 8.d0*&
                      &Absw22(f)*Rew1(f)*Rew2(f) + 4.d0*(Absw12(f) -&
                      &3.d0*Absw22(f))*Rew2(f)**2 - 16.d0*Rew1(f)*&
                      &Rew2(f)**3 + 16.d0*Rew2(f)**4)/(Absw22(f)**2*&
                      &(Absw12(f)**2 + Absw22(f)*(Absw22(f) + 4.d0*&
                      &Rew1(f)*(Rew1(f) - Rew2(f))) - 2.d0*Absw12(f)*&
                      &(Absw22(f) + 2.d0*(Rew1(f) - Rew2(f))*&
                      &Rew2(f)))), f=3,6) /),&
            &al6 = (/ ((-2.d0*(Absw22(f)*(Rew1(f) - 2.d0*Rew2(f)) +&
                      &Rew2(f)*(Absw12(f) + 4.d0*Rew2(f)*(-Rew1(f) +&
                      &Rew2(f)))))/(Absw22(f)**2*(Absw12(f)**2 +&
                      &Absw22(f)*(Absw22(f) + 4.d0*Rew1(f)*(Rew1(f) -&
                      &Rew2(f))) - 2.d0*Absw12(f)*(Absw22(f) + 2.d0*&
                      &(Rew1(f) - Rew2(f))*Rew2(f)))), f=3,6) /)

!     mass anomalous dimension coefficients gamma_1 to gamma_4
         real(dp), dimension(3:6), parameter ::&
            &ga_1 = (/ (2.d0, f=3,6) /),&
            &ga_2 = (/ (101.d0/12.d0 - 5.d0/18.d0*f, f=3,6) /),&
            &ga_3 = (/ (1249.d0/32.d0 - (277.d0/108.d0+5.d0/3.d0*z3)*f&
                       &- 35.d0/648.d0*f*f, f=3,6) /),&
            &ga_4 = (/ (4603055.d0/20736.d0+1060.d0/27.d0*z3-275.d0/&
                       &4.d0*z5 + (11.d0/144.d0*pi**4-91723.d0/3456.d0-&
                       &2137.d0/72.d0*z3+575.d0/36.d0*z5)*f + (2621.d0/&
                       &15552.d0-pi**4/216.d0+25.d0/36.d0*z3)*f*f -&
                       &(83.d0/7776.d0-z3/54.d0)*f**3, f=3,6) /)

!     Auxiliary constants for running quark mass
         real(dp), dimension(3:6), parameter ::&
            &gg_0 = (/ (ga_1(f)/be_1(f), f=3,6) /),&
            &gg_1 = (/ (ga_2(f)/be_1(f)-be_2(f)*ga_1(f)/be_1(f)**2,&
                       &f=3,6) /),&
            &gg_2 = (/ (ga_3(f)/be_1(f)-be_2(f)*ga_2(f)/be_1(f)**2+&
                       &ga_1(f)/be_1(f)**2*(be_2(f)**2/be_1(f)-&
                       &be_3(f)), f=3,6) /),&
            &gg_3 = (/ (ga_4(f)/be_1(f)-be_2(f)*ga_3(f)/be_1(f)**2+&
                       &ga_2(f)/be_1(f)**2*(be_2(f)**2/be_1(f)-&
                       &be_3(f))+ga_1(f)/be_1(f)**2*(be_2(f)*be_3(f)/&
                       &be_1(f)-be_2(f)/be_1(f)*(be_2(f)**2/be_1(f)-&
                       &be_3(f))-be_4(f)), f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &ggg_2 = (/ ((gg_2(f)+gg_1(f)**2)/2.d0, f=3,6) /),&
            &ggg_3 = (/ ((gg_3(f)+1.5d0*gg_1(f)*gg_2(f)+&
                         &gg_1(f)**3/2.d0)/3.d0, f=3,6) /)

         real(dp), dimension(3:6), parameter ::&
            &qm_1 = (/ (ga_1(f)/ga_4(f)-be_ratio_14(f), f=3,6) /),&
            &qm_2 = (/ (ga_2(f)/ga_4(f)-be_ratio_24(f), f=3,6) /),&
            &qm_3 = (/ (ga_3(f)/ga_4(f)-be_ratio_34(f), f=3,6) /)

      end module rge_const
