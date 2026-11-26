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
GEN = cputim.o zeroin.o romb.o splines.o adaptive_matsubara.o spline.o reallocate.o zinv.o \
      c3fft.o ifflim.o ffts.o deter.o dvbes1.o gaunt.o \
      fermi.o notri.o errclr.o errflg.o outerr.o \
      sphbes.o ylm.o stern.o sortag.o readselfenergy.o eigsys.o dfrad.o sorttodiag.o inverss.o

# These routines are tightly coupled, mostly through modules.f90
OBJSMOD =  modules.o structure.o stress_modules.o dmfts.o forces.o mod_mpi.o stress_mpi.o t3j.o rotdef.o stress_rela_operator.o w_atpar.o angle.o charge.o \
           p_project.o csplit.o d5splt.o f7splt.o gtfnam.o \
           harmon.o kdelta.o lomain.o l2mainSubs.o alm_blm_clm.o stress_g2_matrices.o \
	   stress_integrals.o stress_kernels.o stress_rela.o stress_interstitial.o stress_symmetry.o \
	   cmp_mt_density.o cmp_dmft_stress.o cmp_dmft_weights.o frozen_response_io.o l2main.o dmft2.o init_energy_file.o lda_fermi.o \
           latgen.o p3splt.o psplit.o recfil.o recpr.o \
           spin_rotate.o cmp_log_gdloc.o print_pre_atpar.o fomai3.o symmetrize.o

OBJS = $(GEN) $(OBJSMOD)
#.....................................
#  Builds executables and py-modules
#.....................................
all : latpar dmft2 combineud rdVec.so fermif

latpar :
	${MAKE} -C atpar libatpar.a response

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

# Module interfaces introduced by the complete local G2 implementation.
stress_g2_matrices.o: modules.o structure.o w_atpar.o mod_mpi.o
cmp_dmft_stress.o: stress_g2_matrices.o alm_blm_clm.o
l2main.o: stress_g2_matrices.o

# Module interface for the shared local-orbital/projector strain geometry.
l2mainSubs.o alm_blm_clm.o: lomain.o

adaptive_matsubara.o: splines.o
l2main.o cmp_dmft_weights.o: adaptive_matsubara.o
# The self-energy embedding routines and the Fermi module live in l2mainSubs.
l2main.o cmp_dmft_weights.o: l2mainSubs.o
init_energy_file.o l2mainSubs.o dmft2.o: mod_mpi.o

# Scalar-relativistic stress operators share the lifetime of the ATPAR basis.
stress_rela_operator.o: modules.o
w_atpar.o: stress_rela_operator.o
stress_rela.o: w_atpar.o stress_rela_operator.o

# Compare the original radial integration with the default cache on two
# inequivalent atoms, revisiting each with different complex coefficients.
.PHONY: check-rela-operator
check-rela-operator: dmft2
	@set -e; trial_dir=$$(mktemp -d); trap 'rm -rf "$$trial_dir"' EXIT; \
	$(FORT) $(FFLAGS) -I./atpar -o "$$trial_dir/check" tests/test_rela_operator.f90 \
	    $(filter-out dmft2.o,$(OBJS)) $(WLDFLAGS) $(WLIBS) -L./atpar -latpar; \
	export WIEN_STRESS_RELA_DIAG=0 WIEN_STRESS_RELA_OLD_SMALL_SCALE=0.75 \
	    WIEN_STRESS_RELA_SMALL_SURF_SCALE=1 WIEN_STRESS_RELA_ENSMALL_SCALE=1 \
	    WIEN_STRESS_RELA_ENNAT_SCALE=1; \
	for mass in 0 1; do \
	    export WIEN_STRESS_RELA_MASSV_SCALE=$$mass; \
	    EDMFT_STRESS_RELA_REFERENCE=1 "$$trial_dir/check" write "$$trial_dir/reference"; \
	    EDMFT_STRESS_RELA_REFERENCE=0 "$$trial_dir/check" check "$$trial_dir/reference"; \
	done
