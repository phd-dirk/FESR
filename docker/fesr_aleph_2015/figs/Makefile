#
# General purpose Makefile
#

FFLAGS = -Wall -O -L/usr/local/cern/2006/lib

.f.o: ; gfortran $(FFLAGS) -c $*.f
.c.o: ; gcc $(CFLAGS) -c $*.c

OBJECTS1 = main.o
OBJECTS2 = 
OBJECTS = $(OBJECTS1) $(OBJECTS2)

main: $(OBJECTS) Makefile
	gfortran -o main $(OBJECTS) $(FFLAGS)  -lpacklib  -lmathlib -lkernlib

clean: 
	rm *.o main

