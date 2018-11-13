      ! Definition of different sets of s_0's to be used in the analysis
      module s0s_sets
      use num_const
      use params

      implicit none

!     One s_0 at Mtau^2
      real(dp), dimension(0:1), parameter :: s0_mtau = (/ &
      &1.d0, stau /)

!     2 s_0's from 2.9000 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:2), parameter :: s0_2_aleph = (/ &
      &2.d0,&
      &    stau, 3.0000d0 /)

!     3 s_0's from 2.7000 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:3), parameter :: s0_3_aleph = (/ &
      &3.d0,&
      &    stau, 3.0000d0, 2.8000d0 /)

!     4 s_0's from 2.5000 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:4), parameter :: s0_4_aleph = (/ &
      &4.d0,&
      &    stau, 3.0000d0, 2.8000d0, 2.6000d0 /)

!     5 s_0's from 2.3500 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:5), parameter :: s0_5_aleph = (/ &
      &5.d0,&
      &    stau, 3.0000d0, 2.8000d0, 2.6000d0, 2.4000d0 /)

!     6 s_0's from 2.2500 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:6), parameter :: s0_6_aleph = (/ &
      &6.d0,&
      &    stau, 3.0000d0, 2.8000d0, 2.6000d0, 2.4000d0,&
      &2.3000d0 /)

!     7 s_0's from 2.1500 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:7), parameter :: s0_7_aleph = (/ &
      &7.d0,&
      &    stau, 3.0000d0, 2.8000d0, 2.6000d0, 2.4000d0,&
      &2.3000d0, 2.2000d0 /)

!     8 s_0's from 2.0500 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:8), parameter :: s0_8_aleph = (/ &
      &8.d0,&
      &    stau, 3.0000d0, 2.8000d0, 2.6000d0, 2.4000d0,&
      &2.3000d0, 2.2000d0, 2.1000d0 /)

!     9 s_0's from 1.9750 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:9), parameter :: s0_9_aleph = (/ &
      &9.d0,&
      &    stau, 3.0000d0, 2.8000d0, 2.6000d0, 2.4000d0,&
      &2.3000d0, 2.2000d0, 2.1000d0, 2.0000d0 /)

!     43 s_0's from 0.9875 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:43), parameter :: s0_43_aleph = (/ &
      &43.d0,&
      &    stau, 3.000d0, 2.800d0, 2.600d0, 2.400d0, 2.300d0,&
      & 2.200d0, 2.100d0, 2.000d0, 1.950d0, 1.900d0, 1.850d0,&
      & 1.800d0, 1.750d0, 1.700d0, 1.675d0, 1.650d0, 1.625d0,&
      & 1.600d0, 1.575d0, 1.550d0, 1.525d0, 1.500d0, 1.475d0,&
      & 1.450d0, 1.425d0, 1.400d0, 1.375d0, 1.350d0, 1.325d0,&
      & 1.300d0, 1.275d0, 1.250d0, 1.225d0, 1.200d0, 1.175d0,&
      & 1.150d0, 1.125d0, 1.100d0, 1.075d0, 1.050d0, 1.025d0,&
      & 1.000d0 /)

!     78 s_0's from 0.1000 to 3.0875 GeV^2 at upper edge of ALEPH bins
      real(dp), dimension(0:78), parameter :: s0_78_aleph = (/ &
      &78.d0,&
      &    stau, 3.000d0, 2.800d0, 2.600d0, 2.400d0, 2.300d0,&
      & 2.200d0, 2.100d0, 2.000d0, 1.950d0, 1.900d0, 1.850d0,&
      & 1.800d0, 1.750d0, 1.700d0, 1.675d0, 1.650d0, 1.625d0,&
      & 1.600d0, 1.575d0, 1.550d0, 1.525d0, 1.500d0, 1.475d0,&
      & 1.450d0, 1.425d0, 1.400d0, 1.375d0, 1.350d0, 1.325d0,&
      & 1.300d0, 1.275d0, 1.250d0, 1.225d0, 1.200d0, 1.175d0,&
      & 1.150d0, 1.125d0, 1.100d0, 1.075d0, 1.050d0, 1.025d0,&
      & 1.000d0, 0.975d0, 0.950d0, 0.925d0, 0.900d0, 0.875d0,&
      & 0.850d0, 0.825d0, 0.800d0, 0.775d0, 0.750d0, 0.725d0,&
      & 0.700d0, 0.675d0, 0.650d0, 0.625d0, 0.600d0, 0.575d0,&
      & 0.550d0, 0.525d0, 0.500d0, 0.475d0, 0.450d0, 0.425d0,&
      & 0.400d0, 0.375d0, 0.350d0, 0.325d0, 0.300d0, 0.275d0,&
      & 0.250d0, 0.225d0, 0.200d0, 0.175d0, 0.150d0, 0.125d0/)

      end module s0s_sets

