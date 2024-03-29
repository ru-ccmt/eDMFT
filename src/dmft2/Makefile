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
GEN = cputim.o zeroin.o romb.o splines.o spline.o reallocate.o zinv.o \
      c3fft.o ifflim.o ffts.o deter.o dvbes1.o gaunt.o \
      fermi.o notri.o errclr.o errflg.o outerr.o \
      sphbes.o ylm.o stern.o sortag.o readselfenergy.o eigsys.o dfrad.o sorttodiag.o inverss.o

# These routines are tightly coupled, mostly through modules.f90
OBJSMOD =  modules.o structure.o dmfts.o forces.o mod_mpi.o t3j.o rotdef.o w_atpar.o angle.o charge.o \
           p_project.o csplit.o d5splt.o f7splt.o gtfnam.o \
           harmon.o kdelta.o lomain.o l2mainSubs.o alm_blm_clm.o \
	   cmp_mt_density.o cmp_dmft_weights.o l2main.o dmft2.o init_energy_file.o lda_fermi.o \
           latgen.o p3splt.o psplit.o recfil.o recpr.o \
           spin_rotate.o cmp_log_gdloc.o print_pre_atpar.o fomai3.o symmetrize.o

OBJS = $(GEN) $(OBJSMOD)
#.....................................
#  Builds executables and py-modules
#.....................................
all : latpar dmft2 combineud rdVec.so fermif

latpar :
	cd atpar && ${MAKE} libatpar.a

dmft2 : $(OBJS)
	$(FORT) $(LDFLAGS) -o dmft2 $(OBJS) $(WLDFLAGS) $(WLIBS) -L./atpar -latpar

combineud : combineud.f90
	$(F90) -o combineud combineud.f90

fermif : fermi_from_data.o zeroin.o
	$(FORT) $(LDFLAGS) -o $@ fermi_from_data.o zeroin.o

rdVec.so : rdVec.f90
	$(CMP) $(F2PL) -c $? -m rdVec
	mv rdVec.*so rdVec.so

clean:
	rm -f *.o *_tmp_.* *.P .real .complex .sequential .parallel *.mod dmft2 combineud *.so *__genmod.f90 mod_mpi.f90 fermif
	rm -rf *.dSYM
	cd atpar; make clean; cd ../

#...................
# Pattern rules
#...................
%.o : %.f90
	$(FORT) $(FFLAGS) -I./atpar -c $< -o $@

%.o : %.FPP
	$(PF90) $(Fmpi_define)  $(FFLAGS) -c $< -o $@

%.o : %.c
	$(CC) -c $< -o $@

%.so : %.f90
	$(CMP) -c $< -m $*
