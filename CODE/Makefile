ifdef COMPILE
ifeq ($(COMPILE),pathscale)
override define COMPILE
pathf95
endef
TARCH  = -r8 -m64 -march=opteron -mcmodel=medium -ipo -convert little_endian
F95OPT = -intrinsic=PGI -OPT:Ofast
FIXED  = -fixedform
FREE   = -freeform
ifeq ($(ENDIAN),yes)
override define TARCH
-r8 -m64 -march=opteron -mcmodel=medium -convert big_endian
endef
endif
ifeq ($(OMP),no)
TARGET = chymera-ser
MP     =
else
TARGET = chymera-omp
MP     = -mp
endif
CC     = pathcc $(TARCH)
endif
ifeq ($(COMPILE),intel)
override define COMPILE
ifort
endef
TARCH  = -r8 -mcmodel=medium -convert little_endian 
CARCH  = -mcmodel=medium
F95OPT = -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -shared-intel 
FIXED  = -fixed
FREE   = -free
INLINE = -ip
ifeq ($(ENDIAN),yes)
override define TARCH
-r8 -mcmodel=medium -convert big_endian
endef
endif
ifeq ($(OMP),no)
TARGET = chymera-ser
MP     = 
else
TARGET = chymera-omp
MP     = -openmp
endif
CC     = icc $(MP) $(CARCH)
endif
ifeq ($(COMPILE),portland)
override define COMPILE
pgf95
endef
TARCH  = -r8 -tp k8-64
F95OPT = -g -O3
FIXED  = -Mfixed
FREE   = -Mfree
ifeq ($(ENDIAN),yes)
TARCH := $(TARCH) -byteswapio
endif
ifeq ($(OMP),no)
TARGET = chymera-ser
MP     = 
else
TARGET = chymera-omp
MP     = -mp
endif
CC     = pgcc $(MP) $(TARCH)
endif
ifeq ($(COMPILE),gcc)
override define COMPILE
gfortran
endef
TARCH  = -m64 -fconvert=little-endian
F95OPT = -g -O3
FIXED  = -ffixed-form
FREE   = -ffree-form
ifeq ($(ENDIAN),yes)
override define TARCH
-m64 -fconvert=big-endian
endef
endif
ifeq ($(OMP),no)
TARGET = chymera-ser
MP     = 
else 
TARGET = chymera-omp
MP     = -fopenmp
endif
CC = gcc $(TARCH)
endif
else
COMPILE     = gfortran
TARCH  = -m64 -fconvert=big-endian
F95OPT = -g -O3
FIXED  = -ffixed-form
FREE   = -ffree-form
ifeq ($(ENDIAN),yes)
override define TARCH
-m64 -fconvert=little-endian
endef
endif
ifeq ($(OMP),no)
TARGET = chymera-ser
MP     = 
else 
TARGET = chymera-omp
MP     = -fopenmp
endif
CC = gcc $(TARCH) 
endif


OBJ = 3dhyd-main.o fft.o flux.o housekeeping.o io.o misc.o pot3.o \
	rad.o source.o state.o boundary.o blktri.o avisc.o hybrid.o \
	initengtable.o wiggle.o planet.o first_touch.o
FLUIDOBJ  =  fluid_mod.o fluid_start.o fluid_setup.o fluid_interp.o\
             fluid_advance.o fluid_writeout.o fluid_restart.o
OBJINTERP = interp_project.o interp.o
CLINK     = -L$(HOME)/lib/lib -I$(HOME)/lib/include -lgsl -lgslcblas

OBJ         := $(OBJ) $(FLUIDOBJ) $(OBJINTERP)
LFLAGS      := $(TARCH)
CFLAGS      := $(TARCH) $(F95OPT) $(FIXED)
CFLAGS-FREE := $(TARCH) $(F95OPT) $(FREE)

target: $(OBJ) 
	$(COMPILE) $(MP) $(LFLAGS) -o $(TARGET) $(OBJ) $(CLINK)
clean:
	rm -f *.o *.lst chymera-*

#==============================================================================
#                          OBJECTS

3dhyd-main.o: 3dhyd-main.f
	$(COMPILE) $(MP) $(CFLAGS) -c $<

fft.o: fft.f
	$(COMPILE) $(MP) $(CFLAGS) -c $<

flux.o: flux.f
	$(COMPILE) $(MP) $(CFLAGS) $(INLINE)  -c  $<

housekeeping.o: housekeeping.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

io.o: io.f
	$(COMPILE) $(CFLAGS)  -c $<

misc.o: misc.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

pot3.o: pot3.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

rad.o: rad.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

source.o: source.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

state.o: state.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

boundary.o: boundary.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

blktri.o: blktri.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

avisc.o: avisc.f
	$(COMPILE) $(MP) $(CFLAGS)  -c $<

hybrid.o: hybrid.f
	$(COMPILE) $(MP) $(CFLAGS)  -c hybrid.f

initengtable.o: initengtable.f
	$(COMPILE) $(MP) $(CFLAGS)  -c initengtable.f

wiggle.o: wiggle.F
	$(COMPILE) $(MP) $(CFLAGS-FREE)  -c wiggle.F

defineme.o:defineme.f
	$(COMPILE) $(MP) $(CFLAGS) -c defineme.f

interp.o: interp.c
	$(CC) $(COPT) $(CLINK) -c interp.c

planet.o: planet.f
	$(COMPILE) $(MP) $(CFLAGS) -c planet.f

first_touch.o: first_touch.f
	$(COMPILE) $(MP) $(TARCH) $(FIXED) -c first_touch.f
%.o:%.F
	$(COMPILE) $(MP) $(CFLAGS) -c $<
