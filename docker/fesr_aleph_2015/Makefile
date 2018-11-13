#
# Makefile for "fesr" project
# Last change: Matthias Jamin, 12.11.2014
#

FFLAGS = -ffree-form -Wall -O -L/fesr/cernlib/cernlib.web.cern.ch/cernlib/download/2006b_x86_64-slc5-gcc43-opt/lib
# FFLAGS = -ffree-form -Wall -O -L/usr/local/cern/2006/lib
# FFLAGS = -pedantic -Wall -O -L/usr/local/cern/2006/lib

%.o : %.f90 ; gfortran $(FFLAGS) -c $*.f90
# .c.o: ; gcc $(CFLAGS) -c $*.c

OBJECTS1 = num_const.o gauss.o params.o paramsDV.o rge_const.o condens.o \
           weights_D.o weights_rho.o weights_used.o s0s_sets.o s0s_fit.o \
           s0s_out.o cVA0_const.o cVA0_depen.o RPM_const.o opal98_spec.o \
           aleph_vv.o aleph_aa.o aleph_vpa.o aleph_vma.o matrix.o
OBJECTS2 = fesr.o gauleg.o alphas_run.o mq_run.o snleq.o VAcorrD0.o \
           VAcorrD2.o VAcorrD4.o VAcorrD68.o RPMmodel.o incgam.o \
           VAcintD0.o VAcintD2.o VAcintD4.o VAcintD68.o VAdeltas.o \
           PSpheno.o VAmomAleph2014.o VAmomOpal2012.o VAmomth.o dVAmomth.o \
           pipole.o dualvio.o chi2DV.o rtauVpAfit.o
OBJECTS = $(OBJECTS1) $(OBJECTS2)

main: $(OBJECTS) Makefile
	gfortran -o main $(FFLAGS) $(OBJECTS) -lpacklib

clean: 
	rm $(OBJECTS) *.mod # *.out

