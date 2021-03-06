# -*-Makefile-*-
# makefile for the 3-D astrophysical hydrodynamics code.

 #MAKEFILE = Makefile.intel

#===============================================================================
 TARCH      = -r8 -axsse3 -mcmodel=medium -convert big_endian

 OMPF95     = ifort

 F90FLAGS   = -O3 -pad -align all -ipo -no-prec-div -nocheck -shared-intel -openmp -fixed

 LIBS       = $(F90FLAGS) $(TARCH)

# C compilers
 CC         = icc

 .f.o:
	$(OMPF95) $(F90FLAGS) -c $*.f

 OBJ = fft.o dcomp.o blktri.o torqueout.o boundary.o initengtable.o

#===============================================================================
#                        EXECUTABLES

decompose: dummy
	make target "FC=$(OMPF95)" "CC=$(CC)" \
 "CFLAGS=$(F90FLAGS) $(TARCH)" "TARGET = decompose"

#-------------------------------------------------------------------------------
target: $(OBJ) $(PERFOBJ) $(DBGOBJ)
	$(FC) $(LIBS) -o $(TARGET) $(OBJ)

dummy:

clean:
	rm -f *.o decompose


#===============================================================================
#                          OBJECTS

fft.o: fft.f
	$(FC) $(CFLAGS) -c $<

dcomp.o: dcomp.f
	$(FC) $(CFLAGS) -c dcomp.f

blktri.o: blktri.f
	$(FC) $(CFLAGS) -c $<

torqueout.o: torqueout.F
	$(FC) $(CFLAGS) -c $<

boundary.o: boundary.f
	$(FC) $(CFLAGS) -c $<
initengtable.o: initengtable.f
	$(FC) $(CFLAGS) -c $<
