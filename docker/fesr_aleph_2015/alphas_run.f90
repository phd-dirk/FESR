
!     Routines for running of alpha_s in the complex plane
!     Last change: Matthias Jamin, 2.3.2009

!-----------------------------------------------------------------------

!     Calculates a(q^2) from integrating the RG-equation
!     in the complex q^2 plane from a given a(p^2) at p^2

      function zarg(nf,q2,p2,ap)
      use num_const

      implicit none
      integer,     intent(in) :: nf
      real(dp),    intent(in) :: p2, ap
      complex(dc), intent(in) :: q2
      complex(dc)             :: zarg

      integer     :: info, f
      real(dp)    :: xx(2), ff(2), ww(10)
      complex(dc) :: zaint4, zaint5, zaex, aprox, arem
      common /ACOM/ arem, f
      external asub
       
      f = nf

!     Select 4-loop or 5-loop running
!      arem  = zaint4(f,cmplx(ap,0,dc))-0.5d0*log(q2/p2)
     arem  = zaint5(f,cmplx(ap,0,dc))-0.5d0*log(q2/p2)
       
      aprox = zaex(nf,q2)
      xx(1) = dreal(aprox)
      xx(2) = aimag(aprox)
       
!     My version of the routine to solve non-linear equations
      call my_snleq(2,xx,ff,1.D-15,1.D-15,300,0,info,asub,ww)
!     CERNlib version of the routine to solve non-linear equations
!     call   dsnleq(2,xx,ff,1.D-15,1.D-15,300,0,info,asub,ww)
!     Write (6,*) cmplx(F(1),F(2),dc), info
       
      zarg = cmplx(xx(1),xx(2),dc)
       
      return
      end function zarg
!-----------------------------------------------------------------------

!     Subroutine to provide the set of non-linear equations
!     for calculating alpha_s(q^2) in the complex plain

      subroutine asub(n,xx,ff,k)
      use num_const

      implicit none
      integer,  intent(in)  :: n, k
      real(dp), intent(in)  :: xx(n)
      real(dp), intent(out) :: ff(n)

      integer     :: f
      complex(dc) :: aeq, arem, zaint4, zaint5
      common /ACOM/ arem, f

!     Select 4-loop or 5-loop running
!      aeq = zaint4(f,cmplx(xx(1),xx(2),dc))-arem
     aeq = zaint5(f,cmplx(xx(1),xx(2),dc))-arem

      if (k == 1) then
         ff(1) = dreal(aeq)
         else if (k.eq.2) then
         ff(2) = aimag(aeq)
         else
         write(*,*) "Error in ASUB!"
      end if

      return
      end subroutine asub
!-----------------------------------------------------------------------

!     4-loop integral of the RGE for alpha_s(q)/pi   (complex version)

      function zaint4(nf,x)
      use num_const
      use rge_const

      implicit none
      integer,     intent(in) :: nf
      complex(dc), intent(in) :: x
      complex(dc)             :: zaint4

      complex(dc) :: IsqrtD
       
      f = nf
      IsqrtD = ii*sqrt(4*qw_2(f)-qw_1(f)**2)

      zaint4 = ((qw_2(f)-qw_1(f)*(qw_1(f)-qw_3(f)))/2/qw_2(f)**2*&
              &log(x*x/(x*x+qw_1(f)*x+qw_2(f)))-(qw_1(f)-qw_3(f))/&
              &qw_2(f)/x-(qw_1(f)*qw_2(f)+(qw_1(f)-qw_3(f))*(2*qw_2(f)-&
              &qw_1(f)**2))/2/qw_2(f)**2/IsqrtD*log((IsqrtD-qw_1(f)-2*&
              &x)/(IsqrtD+qw_1(f)+2*x))+1/qw_3(f)/x-log((qw_3(f)+x)/x)&
              &/qw_3(f)**2)/(qw_3(f)*(qw_1(f)-qw_3(f))-qw_2(f))/be_4(f)

      return
      end function zaint4
!-----------------------------------------------------------------------

!     5-loop integral of the RGE for alpha_s(q)/pi   (complex version)

      function zaint5(nf,x)
      use num_const
      use rge_const

      implicit none
      integer,     intent(in) :: nf
      complex(dc), intent(in) :: x
      complex(dc)             :: zaint5

      complex(dc) :: W1, W2

      f = nf
      W1 = i1*sqrt(Absw12(f) - Rew1(f)**2)
      W2 = i1*sqrt(Absw22(f) - Rew2(f)**2)

      zaint5 = ( al2(f)*log(x) - al1(f)/x +&
     &  al4(f)/2.d0*log(x*x - 2.d0*Rew1(f)*x + Absw12(f)) +&
     &  (al3(f) + al4(f)*Rew1(f))/W1*atan((x - Rew1(f))/W1) +&
     &  al6(f)/2.d0*log(x*x - 2.d0*Rew2(f)*x + Absw22(f)) +&
     &  (al5(f) + al6(f)*Rew2(f))/W2*atan((x - Rew2(f))/W2) )/be_5(f)

      return
      end function zaint5
!-----------------------------------------------------------------------

!     Expanded version of 4-loop expression for a(q^2)   (complex)

      function zaex(nf,q2)
      use num_const
      use rge_const
      use params

      implicit none
      integer,     intent(in) :: nf
      complex(dc), intent(in) :: q2
      complex(dc)             :: zaex

      complex(dc) :: L, lL

      f = nf
      L  = log(q2/lambda(f)**2)
      lL = log(L)

      zaex = 2/be_1(f)/L*(1-2*be_ratio_21(f)*lL/be_1(f)/L+&
            &4/(be_1(f)*L)**2*(be_ratio_21(f)**2*(lL**2-lL-1)+&
            &be_ratio_31(f))+8/(be_1(f)*L)**3*(be_ratio_21(f)**3*&
            &((-lL)**3+2.5d0*lL**2+2*lL-0.5d0)-3*be_ratio_21(f)*&
            &be_ratio_31(f)*lL+be_ratio_41(f)/2))

      return
      end function zaex
!-----------------------------------------------------------------------
