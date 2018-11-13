!     Definition of different sets of s_0's to be used in the analysis
      module s0s_out
      use num_const
      use weights_used

         implicit none
         integer, parameter :: Nsum = 78 ! Total number of s_0's

         integer,  dimension(0:Nmom) :: Ns0s    ! # s_0's per moment
         real(dp), dimension(Nmom,0:80) :: s0s ! Sets of used s_0's

      contains

!     Initialise set of s_0's

         subroutine init_s0s_out(Nsum,Ns0s,s0s)
         use weights_used
         use s0s_sets

         implicit none
         integer :: Nsum
         integer,  dimension(0:Nmom) :: Ns0s
         real(dp), dimension(Nmom,0:80) :: s0s

         integer  :: n, sun

!      Used sets of s_0's
!      Make sure that the number of definitions = Nmom !
!           s0s(1,0:1) = s0_mtau
!           s0s(1,0:5) = s0_5_aleph
            s0s(1,0:78) = s0_78_aleph
!           s0s(2,0:1) = s0_mtau
!           s0s(3,0:1) = s0_mtau
!           s0s(4,0:1) = s0_mtau
!           s0s(5,0:1) = s0_mtau

!      Set numbers of s_0's
         sun = 0
         do n = 1, Nmom
            Ns0s(n) = int(s0s(n,0))
            sun = sun + Ns0s(n)
         end do

         Ns0s(0) = sun

         if (Ns0s(0) /= Nsum) then
            write(*,*) "Inconsistent total number of s_0's!"
            stop
         end if

      end subroutine init_s0s_out

      end module s0s_out
