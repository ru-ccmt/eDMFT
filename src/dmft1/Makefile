include ../Makefile.in

FORT = $(PF90)
LDFLAGS= $(WFOPT)
FFLAGS  = $(WFOPT)
###############################################################################

#..............................................................................
#
#  Object files 
#
# These routines are completely general and decoupled from particular code. Can be reused in other parts of the code
#..............................................................................
GEN = cputim.o splines.o stern.o

# These routines are tightly coupled, mostly through modules.f90
OBJSMOD = mod_mpi.o error_handle.o structure.o gtfnam2.o modules.o w_atpar.o angle.o atpar.o dvbes1.o \
          romb.o harmon.o inverss.o spline.o p_project.o cmpdmftrans.o printgloc.o l2mainSubs.o l2main.o eigsys.o \
          dmftmain.o  rotate.o sphbes.o symop.o ylm.o mknam.o \
          sym.o lomain.o cmpxqtl.o testxqtldiff.o readselfenergy.o \
	  zinv.o overlap1.o spin_rotate.o cmp_dm.o symmetrize.o

# Python f2py routines
pylibs = trafo.so trafoso.so occupr.so occupi.so rdU.so

OBJS = $(GEN) $(OBJSMOD)
#..............................................................................

all : latpar dmft $(pylibs)

dmft : $(OBJS)
	$(FORT) $(LDFLAGS) -o dmft $(OBJS) $(WLDFLAGS) $(WLIBS) -L./atpar -latpar

latpar :
	cd atpar && ${MAKE} libatpar.a

trafo.so : trafo.f90
	$(CMP) -c trafo.f90 -m trafo
	mv trafo.*so trafo.so

trafoso.so : trafoso.f90
	$(CMP) -c trafoso.f90 -m trafoso
	mv trafoso.*so trafoso.so

occupr.so : occupr.f90
	$(CMP) -c occupr.f90 -m occupr
	mv occupr.*so occupr.so

occupi.so : occupi.f90
	$(CMP) -c occupi.f90 -m occupi
	mv occupi.*so occupi.so

rdU.so : rdU.f90
	$(CMP) -c rdU.f90 -m rdU
	mv rdU.*so rdU.so

#..............................................................................
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod *.so dmft mod_mpi.f90 atpar/libatpar.a
	rm -rf *.dSYM
	cd atpar; make clean
#..............................................................................
#
#  define inference rules to generate object files from source files
#
#...................
# Pattern rules
#...................
%.o : %.f90
	$(FORT) $(FFLAGS) -c $< -o $@

%.o : %.F90
	$(FORT) $(FFLAGS) -c $< -o $@

%.o : %.FPP
#	$(preproc) -P  $(Fmpi_define) $< $*.f90
#	$(PF90) $(FFLAGS) -c $*.f90
	$(PF90) $(Fmpi_define) $(FFLAGS) -c $< -o $@
#	rm $*.f90

%.o : %.c
	$(CC) -c $< -o $@

%.so : %.f90
	$(CMP) -c $< -m $*


