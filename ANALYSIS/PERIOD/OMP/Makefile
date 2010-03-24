#   $Id: Makefile,v 1.11 2005/02/01 20:24:11 scottm Exp $
#
#  Copyright (c) 1988-2005, Research Systems Inc.  All rights reserved.
#  This software includes information which is proprietary to and a
#  trade secret of Research Systems, Inc.  It is not to be disclosed
#  to anyone outside of this organization. Reproduction by any means
#  whatsoever is  prohibited without express written permission.
#
#  Platform independant makefile for the IDL Call_External() examples.
#
SHELL=/bin/sh
CFLAGS=
FFLAGS=

SHL_LIBS=  global_coefs.$(SHL_EXT) 

# The following is the default entry point. This section will determine 
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

# Note: For Sun Forte versions prior to version 7, change the
#       F_LD_POST below to F_LD_POST= -lF77 -lm -lc

all : 
	@echo "OS type detected: "`uname`
	@case `uname` in \
	   "SunOS") if [  -d /proc ]; then \
			make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"F77=f77" \
			"C_LD=cc" \
			"F_LD=f77" \
			"CFLAGS=-K pic -DSPARC -G" \
			"C_LD_FLAGS=-G -DSOLARIS" \
			"FFLAGS=-pic -G" \
			"FOBJS=ftn_only_sun.o ftn_strarr_sun.o" \
			"F_LD_FLAGS=-G" \
			"F_LD_POST= -lfsu -lm -lc"; \
	   	    fi  \
		;; \
	   "AIX") make libs \
			"SHL_EXT=a" \
			"CC=cc" \
			"F77=f77" \
			"C_LD=cc" \
			"F_LD=f77" \
			"C_LD_FLAGS=-bM:SRE -bnoentry" \
			"F_LD_FLAGS=-bM:SRE -bnoentry" \
			"EXAMPLE_C2F_EXPORT=-bE:example_c2f.exp" \
			"EXAMPLE1_EXPORT=-bE:example_c2f.exp" \
			"SIMPLEC2F_EXPORT=-bE:simple_c2f.exp" \
			"FORTRAN_EXPORT=-bE:ftn_only.exp" \
			"FOBJS=ftn_only_ibm.o" \
			"F_LD_POST= -lxlf -lxlf90";; \
	   "HP-UX") make libs \
			"SHL_EXT=sl" \
			"CC=cc" \
			"F77=f77" \
			"C_LD=ld" \
			"F_LD=ld" \
			"CFLAGS=+z" \
			"C_LD_FLAGS=-b" \
			"F_LD_FLAGS=-b" \
			"FFLAGS=+z" \
			"FOBJS=ftn_only_hp.o ftn_strarr_hp.o" ;; \
	   "IRIX" ) make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"F77=@echo \"NOT Implemented for example: f77 \"" \
			"C_LD=ld" \
			"F_LD=@echo \"NOT Implemented for example: f77 \"" \
			"CFLAGS=-KPIC" \
			"C_LD_FLAGS=-shared";; \
	   "IRIX64" ) make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"F77=@echo \"NOT Implemented for example: f77 \"" \
			"C_LD=ld" \
			"F_LD=@echo \"NOT Implemented for example: f77 \"" \
			"CFLAGS=-KPIC -ansi -32" \
			"C_LD_FLAGS=-shared -32";; \
	   "OSF1" )  make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"F77=f77" \
			"C_LD=ld" \
			"F_LD=ld" \
			"CFLAGS=-DOSF1" \
			"C_LD_FLAGS=-S -shared" \
			"F_LD_FLAGS=-S -shared" \
			"FOBJS=ftn_only_osf1.o" \
			"F_LD_POST= -lUfor -lfor -lFutil -lm -lots -lc";; \
	   "Linux" ) make libs \
			"SHL_EXT=so" \
			"CC=gcc" \
			"F77=pgf90" \
			"C_LD=ld" \
			"F_LD=pgf90" \
			"CFLAGS=-m64 -fPIC -D_REENTRANT" \
			"FFLAGS=-tp k8-64 -Mfree -fPIC -C"\
			"F_LD_FLAGS=-shared"\
			"C_LD_FLAGS=-shared";; \
	   *) echo "This system is not supported" ;; \
       esac

libs : $(SHL_LIBS)
	@echo ""
	@echo "To test the call external examples, start IDL and execute one"
	@echo "of the following IDL procedures: simple_c2f or ftn_only"

# Fortran 

global_coefs.$(SHL_EXT) : globalinclude.o globalwrapper.o
	$(F_LD) $(F_LD_FLAGS) -o global_coefs.$(SHL_EXT) \
		globalinclude.o globalwrapper.o $(EXAMPLE_C2F_EXPORT) $(F_LD_POST)

example1.$(SHL_EXT) : example1.o
	$(F_LD) $(F_LD_FLAGS) -o example1.$(SHL_EXT) \
		example1.o $(EXAMPLE1_EXPORT) $(F_LD_POST)

simple_c2f.$(SHL_EXT) : simple_c2f.o simple_c2f1.o
	$(F_LD) $(F_LD_FLAGS) -o simple_c2f.$(SHL_EXT) \
		simple_c2f.o simple_c2f1.o $(SIMPLEC2F_EXPORT) $(F_LD_POST)


ftn_only.$(SHL_EXT) : $(FOBJS)
	$(F_LD) $(F_LD_FLAGS) -o ftn_only.$(SHL_EXT) $(FOBJS) \
		$(FORTRAN_EXPORT) $(F_LD_POST)

# Cleanup

tidy :
	rm -f *.o

clean :
	rm -f *.o *.so *.sl *.a


#  Rules for building the object files.

.SUFFIXES: .c .o .f

.c.o :
	$(CC) $(CFLAGS) -c $*.c

.f.o :
	$(F77) $(FFLAGS)  -c $*.f


