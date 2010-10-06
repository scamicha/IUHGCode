# -*-Makefile-*-
#===============================================================================
 TARCH      = -r8 -axsse3 -mcmodel=medium #-convert big_endian

 OMPF95     = ifort

 F90FLAGS   = -O3 -pad -align all -ipo -no-prec-div -nocheck -shared-intel -openmp

 LIBS       = $(F90FLAGS) $(TARCH)

 .f.o:
	$(OMPF95) $(F90FLAGS) -c $*.f

 OBJ = qplot_eos.o initengtable.o

#===============================================================================
#                        EXECUTABLES

decompose: dummy
	make target "FC=$(OMPF95)" \
 "CFLAGS=$(F90FLAGS) $(TARCH)" "TARGET = qeos"

#-------------------------------------------------------------------------------
target: $(OBJ) $(PERFOBJ) $(DBGOBJ)
	$(FC) $(LIBS) -o $(TARGET) $(OBJ)

dummy:

clean:
	rm -f *.o qeos


#===============================================================================
#                          OBJECTS

qplot_eos.o: qplot_eos.f90
	$(FC) $(CFLAGS) -c $<
initengtable.o: initengtable.f90
	$(FC) $(CFLAGS) -c $<
