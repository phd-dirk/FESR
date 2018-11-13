      ! Definition of the purely perturbative coefficients
      ! of the vector and axialvector correlators
      module cVA0_const
         use num_const

         implicit none
!     Perturbative coefficients of the non-singlet vector correlator
         real(dp), dimension(0:12,1:12,3:6) :: cV0
         data cV0(0,1,:) /4*1.d0/
         data cV0(1,1,:) /4*1.d0/
         data cV0(2,1,:) /1.63982120489698476474d0,&
                         &1.52452580700338095500d0,&
                         &1.40923040910977714527d0,&
                         &1.29393501121617333554d0/
         data cV0(3,1,:) /6.37101448310094071138d0,&
                         &2.75861606806868110241d0,&
                        &-0.68136860573140256255d0,&
                        &-3.94893953829931028351d0/
         data cV0(4,1,:) /49.07570000294798513221d0,&
!        data cV0(4,1,:) /25.d0,& ! ALEPH report 2005
                         &27.38880009908030194462d0,&
                          &9.21017634921810422202d0,&
                         &-5.52072794676472243940d0/
         data cV0(5,1,3) /283.d0/ ! Guesstimate
!        data cV0(5,1,3) /  0.d0/ ! ALEPH report 2005
!        data cV0(5,1,3) /378.d0/ ! Davier et al. 2008

!        data cV0(6:12,1,3) / 2.91194d3,& ! Davier et al. 2008
!     Higher order coefficients from renormalon model only for n_f = 3
         data cV0(6:12,1,3) / 3.275433197407019d3,&
                             &1.875837290013642d4,&
                             &3.884423798234073d5,&
                             &9.191213065815193d5,&
                             &8.367008582855809d7,&
                            &-5.194684844789161d8,&
                             &3.378627108743168d10/

      end module cVA0_const
