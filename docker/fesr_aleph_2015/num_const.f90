      ! Definition of some numerical constants
      module num_const

         implicit none
         integer, parameter :: dp = kind(1.d0),&
                              &dc = kind((0.d0,1.d0))
      !
         complex(dc), parameter :: i0 = cmplx(0.d0,0.d0,dc),&
                                  &i1 = cmplx(1.d0,0.d0,dc),&
                                  &ii = cmplx(0.d0,1.d0,dc)
      !
         real(dp), parameter    :: pi = 2*asin(1.d0),&
                                  &eg = 0.577215664901532861d0,&
                                  &z3 = 1.202056903159594285d0,&
                                  &z4 = pi**4/90,&
                                  &z5 = 1.036927755143369926d0,&
                                  &z7 = 1.008349277381922827d0
      !
      end module num_const
