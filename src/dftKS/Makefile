include ../Makefile.in

FORT = $(PF90)
LFLIBS = $(WLDFLAGS) $(WLIBS) -L./lattice -L./atpar -llattice -latpar
FFLAGS = $(WFOPT) -I./../lattice -I./../atpar
LDFLAGS = $(WFOPT)

DESTDIR = ./

# We create separate executable for real and complex version.
# The object files (and fortran modules) for real and complex version are compiled into real/  and complex/ subdirectories, respectively.
# Two libraires atpar and lattice are created in the subdirectories with the same name.
#..............................................................................
# common routines, which do not depend on real/complex version
R_COMM = real/ylm.o real/sphbes.o real/dvbes1.o real/gaunt1.o real/vectf.o real/cputim.o real/outerr.o \
	real/find_nloat.o real/cbcomb.o real/lmsort.o real/readFirstPart_of_in1.o real/readSecondPart_of_in1.o  \
	real/readThirdPart_of_in1.o real/readKlist.o
# these routines depend on the type and need preprocessor
R_OBJS = real/dsyr2m.o real/stern.o real/modules.o real/mpi.o real/gtfnam.o  real/errflg.o real/errclr.o real/hamilt.o real/hns.o real/horb.o \
	real/read_and_initialize.o real/lapw1.o real/setwar.o real/tapewf.o real/warpin.o real/wfpnt.o real/wfpnt1.o \
	real/diagonalize.o real/atpar.o real/generate_reciprocal.o real/forfhs.o real/lopw.o real/nn.o real/prtkpt.o \
	real/linEnergySearch.o real/make_albl.o real/lohns.o real/readorbital.o real/read_def_file.o
#..............................................................................
# common routines, which do not depend on real/complex version
C_COMM = cmplx/ylm.o cmplx/sphbes.o cmplx/dvbes1.o cmplx/gaunt1.o cmplx/vectf.o cmplx/cputim.o cmplx/outerr.o \
	cmplx/find_nloat.o cmplx/cbcomb.o cmplx/lmsort.o cmplx/readFirstPart_of_in1.o cmplx/readSecondPart_of_in1.o \
	cmplx/readThirdPart_of_in1.o cmplx/readKlist.o
# these routines depend on the type and need preprocessor
C_OBJS = cmplx/zher2m.o cmplx/stern.o cmplx/modules.o cmplx/mpi.o cmplx/gtfnam.o cmplx/errflg.o cmplx/errclr.o cmplx/hamilt.o cmplx/hns.o cmplx/horb.o \
	cmplx/read_and_initialize.o cmplx/lapw1.o cmplx/setwar.o cmplx/tapewf.o cmplx/warpin.o cmplx/wfpnt.o cmplx/wfpnt1.o \
	cmplx/diagonalize.o cmplx/atpar.o cmplx/generate_reciprocal.o cmplx/forfhs.o cmplx/lopw.o cmplx/nn.o cmplx/prtkpt.o \
	cmplx/linEnergySearch.o cmplx/make_albl.o cmplx/lohns.o cmplx/readorbital.o cmplx/read_def_file.o 

all : latpar llattice lapw1 lapw1c

lapw1 : $(R_COMM) $(R_OBJS) 
	$(FORT) -I./real $(LDFLAGS) -o $@ $(R_COMM) $(R_OBJS) $(LFLIBS)

lapw1c : $(C_COMM) $(C_OBJS) 
	$(FORT) -I./cmplx $(LDFLAGS) -o $@ $(C_COMM) $(C_OBJS) $(LFLIBS)

latpar :
	cd atpar && ${MAKE} libatpar.a
llattice:
	cd lattice && ${MAKE} liblattice.a

clean:
	rm -f real/* cmplx/*
	rm -f lapw1 lapw1c
	rm -f *.mod
	cd atpar; make clean; cd ../
	cd lattice; make clean; cd ../
#..............................................................................
#
#  define inference rules to generate object files from source files
#  (From .FP files the appropriatly typed subroutine (REAL/COMPLEX) has to be chosen.)
#
real/%.o : %.c
	cd real; $(CC) -c ../$< -o $*.o
cmplx/%.o : %.c
	cd cmplx; $(CC) -c ../$< -o $*.o 

real/%.o : %.f90
	cd real; $(FORT) $(FFLAGS) -c ../$< -o $*.o

cmplx/%.o : %.f90
	cd cmplx; $(FORT) $(FFLAGS) -c ../$< -o $*.o

real/%.o : %.FPP
	cd real; $(FORT) -D_REAL $(Fmpi_define) $(FFLAGS) -c ../$< -o $*.o

cmplx/%.o : %.FPP
	cd cmplx; $(FORT) $(Fmpi_define)  $(FFLAGS)  -c ../$< -o $*.o

