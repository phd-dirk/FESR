
!     Routines for running of quark masses in the complex plane
!     Last change: Matthias Jamin, 3.3.2009

!-----------------------------------------------------------------------

!     Calculates ratio m(q^2)/m(p^2) from integrating the RG-equation
!     in the complex q^2 plane from a given a(p^2) at p^2

      function zmrg(nf,q2,p2,atau)
      use num_const
      use params

      implicit none
      integer,     intent(in) :: nf
      real(dc),    intent(in) :: p2, atau
      complex(dc), intent(in) :: q2
      complex(dc)             :: zmrg

      complex(dc) :: ap, aq, zarg, zmint

      ap = zarg(nf,cmplx(p2,0,dc),stau,atau)
      aq = zarg(nf,q2,stau,atau)

      zmrg = exp(zmint(nf,aq)-zmint(nf,ap))

      return
      end
!-----------------------------------------------------------------------

!     Integral of the RGE for m_q(q)   (complex version)

      function zmint(nf,x)
      use num_const
      use rge_const
      use params

      implicit none
      integer,     intent(in) :: nf
      complex(dc), intent(in) :: x
      complex(dc)             :: zmint

      complex(dc) :: J0, J1, J2, J3, IsqrtD

      f = nf
      IsqrtD = ii*sqrt(4.d0*qw_2(f)-qw_1(f)**2)

      J0 = 1.d0/IsqrtD*log((IsqrtD-qw_1(f)-2.d0*x)/&
                          &(IsqrtD+qw_1(f)+2.d0*x))
      J1 = -log((x+qw_3(f))/x)/qw_3(f)
      J2 = (log(x*x/(x*x+qw_1(f)*x+qw_2(f)))-qw_1(f)*J0)/2.d0/qw_2(f)
      J3 = (J0-J1+(qw_1(f)-qw_3(f))*J2)/(qw_3(f)*(qw_1(f)-qw_3(f))-&
          &qw_2(f))

      zmint = ga_4(f)/be_4(f)*(log(x)+qm_3(f)*J1+(qm_2(f)-qm_3(f)*&
             &qw_1(f))*J2+(qm_1(f)-qm_2(f)*qw_3(f)+qm_3(f)*&
             &(qw_1(f)*qw_3(f)-qw_2(f)))*J3)

      return
      end
!-----------------------------------------------------------------------

!     Expanded version of 4-loop expression for m(q^2)   (complex)

      function zmex(nf,q2,p2)
      use num_const
      use rge_const

      implicit none
      integer,     intent(in) :: nf
      real(dc),    intent(in) :: p2
      complex(dc), intent(in) :: q2
      complex(dc)             :: zmex

      complex(dc) :: ap, aq, zaex

      f = nf
      ap = zaex(f,cmplx(p2,0,dc))
      aq = zaex(f,q2)

      zmex = (aq/ap)**gg_0(f)*(1.d0+gg_1(f)*(aq-ap)+ggg_2(f)*(aq**2-&
            &ap**2)-gg_1(f)**2*ap*(aq-ap)+ggg_3(f)*(aq**3-ap**3)-&
            &gg_1(f)*ggg_2(f)*ap*(aq**2+aq*ap-2.d0*ap**2)+&
            &gg_1(f)**3*ap**2*(aq-ap))

      return
      end
!-----------------------------------------------------------------------
