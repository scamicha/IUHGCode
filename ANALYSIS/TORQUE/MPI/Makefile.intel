#-*-Makefile-*-
# makefile for the 3-D astrophysical hydrodynamics code.

 #MAKEFILE = Makefile.intel

#===============================================================================
 TARCH      = -r8 -xHost -mcmodel=medium -convert big_endian

 OMPF95     = mpif90

 F90FLAGS   = -O3 -pad -align all -nocheck -ipo -no-prec-div -shared-intel -openmp -fixed

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

fft.o: fft_mpi.f
	$(FC) $(CFLAGS) -c $<

dcomp.o: dcomp_mpi.f
	$(FC) $(CFLAGS) -c dcomp_mpi.f

blktri.o: blktri_mpi.f
	$(FC) $(CFLAGS) -c $<

torqueout.o: torqueout_mpi.F
	$(FC) $(CFLAGS) -c $<

boundary.o: boundary_mpi.f
	$(FC) $(CFLAGS) -c $<

initengtable.o: initengtable_mpi.f
	$(FC) $(CFLAGS) -c $<
