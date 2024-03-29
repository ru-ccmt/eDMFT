###############################################################################
#                                                                             #
#  make           ... generate executable for the REAL-version                #
#  make complex   ... generate executable for the COMPLEX-version             #
#  make clean     ... delete unnecessary files                                #
#                                                                             #
###############################################################################
#
# FC ........... compiler name
# FOPT ......... optimization flags
# FGEN ......... code generation flags (flags not related to optimization)
# LDFLAGS ...... linker flags
# R_LIBS ....... libraries needed to build the REAL executable
# C_LIBS ....... libraries needed to build the COMPLEX executable
# DESTDIR ...... dir. where the executable should go 
# R_EXECNAME ... name of the resulting REAL executable (without prefixed path!)
# C_EXECNAME ... name of the resulting COMPLEX executable ( -"- )
#
#
include ../../Makefile.in
SHELL = /bin/sh
FC = $(F90) 
CC = cc
FOPT = $(WFOPT) # -ffree-line-length-none -ffree-form -O3 # -FR -DFFTW3 -no-prec-div -DFFTW3 -I/opt/intel/mkl/lib/intel64//include
FGEN = 
LDFLAGS = $(WLIBS) #-framework Accelerate -lfftw3_omp -lfftw3 -lm -L/usr/local/lib #$(FOPT) -L$(MKLROOT)/lib/$(MKL_TARGET_ARCH) -mkl -L/usr/lib64 -pthread
R_LIBS = 
C_LIBS = 
DESTDIR	= ./
R_EXECNAME = optic
C_EXECNAME = opticc

###############################################################################

TYPE='REAL'
TYPE_COMMENT='!_REAL'
FFLAGS = $(FGEN) $(FOPT)
R_EXEC = $(DESTDIR)$(R_EXECNAME)
C_EXEC = $(DESTDIR)$(C_EXECNAME)

#..............................................................................
#
#  Object files common to both REAL and COMPLEX type subroutines
#
OBJS  = abc.o angle.o atpar_op.o  cputim.o  dvbes1.o errclr.o errflg.o gtfnam.o \
        harmon.o  latgen.o  lomain_op.o  mmatcv.o  mmatrix.o  momradintc.o  \
        opmain.o  outerr.o  outmat.o outmatABZ.o outwin.o  planew.o \
        radint.o  rint13.o  rintuout.o rotate.o  \
        sph-UP.o  sph-UPcor.o sphbes.o  symop.o  sym.o  ylm.o

PREOBJS = modules.o
#..............................................................................
#
#  Object files for REAL and COMPLEX versions
#
R_OBJS = $(PREOBJS) $(OBJS)
C_OBJS = $(PREOBJS) $(OBJS)

#..............................................................................
#
#  Some source files contain both REAL and COMPLEX versions. These files
#  have the file extension '.frc' and must be treated differently.
#
.SUFFIXES:	.frc
.SUFFIXES:	.frc90

#..............................................................................
#
#  Build executable (either REAL or COMPLEX versions)
#

all : real complex

real: rm_c_files
	$(MAKE) TYPE='REAL' TYPE_COMMENT='!_REAL' $(R_EXEC)
	rm param.inc

complex: rm_r_files
	$(MAKE) TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' $(C_EXEC)
	rm param.inc

#..............................................................................

$(R_EXEC): $(R_OBJS)
	$(FC) -o $(R_EXEC) $(R_OBJS) $(LDFLAGS) $(R_LIBS)

$(C_EXEC): $(C_OBJS)
	$(FC) -o $(C_EXEC) $(C_OBJS) $(LDFLAGS) $(C_LIBS)

rm_c_files:
	if [ -f .complex ]; then \
	   $(MAKE) clean; \
	fi
	touch .real
	cp -p param.inc_r param.inc

rm_r_files:
	if [ -f .real ]; then \
	   $(MAKE) clean; \
	fi
	touch .complex
	cp -p param.inc_c param.inc

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(R_OBJS):		param.inc modules.frc Makefile
$(C_OBJS):		param.inc modules.frc Makefile

#..............................................................................
#
#  some subroutines require special optimization options on HP-UX 8.05 f77
#
#hns.o:	hns.frc lohns.f rotate.f ylm.f param.inc
#	@sed -e 's/^$(TYPE_COMMENT)//'   hns.frc >_tmp_.f
#	@echo "hns.frc: $(TYPE) version extracted"
#	$(FC) +OP4 -w -WP,-unroll=6,-unroll2=24,-minvector=6 \
#		-WP,-inff='rotate.f',-inff='ylm.f',-inff='lohns.f' -c _tmp_.f
#	mv _tmp_.o hns.o
#
#hamilt.o:	hamilt.frc dvbes1.f sphbes.f rotate.f ustphx.frc \
#		warpin.f ylm.f param.inc
#	$(FC) +OP4c -WP,-unroll=6,-unroll2=24,-minvector=6 -c sphbes.f
#	mv sphbes.P _tmp1_.f
#	@sed -e 's/^$(TYPE_COMMENT)//'   ustphx.frc >_tmp2_.f
#	@echo "ustphx.frc: $(TYPE) version extracted"
#	@sed -e 's/^$(TYPE_COMMENT)//'   hamilt.frc >_tmp_.f
#	@echo "hamilt.frc: $(TYPE) version extracted"
#	$(FC) +OP4 -w -WP,-unroll=6,-unroll2=24,-minvector=6 \
#		-WP,-inff='dvbes1.f',-inff='_tmp1_.f',-inff='rotate.f' \
#		-WP,-inff='ylm.f',-inff='_tmp2_.f',-inff='warpin.f' -c _tmp_.f
#	mv _tmp_.o hamilt.o
#	rm _tmp1_.f _tmp2_.f
#
#sphbes.o:	sphbes.f param.inc
#	$(FC) +OP4 -w -WP,-unroll=6,-unroll2=24,-minvector=6 -c sphbes.f


#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	-rm -f *.o _tmp_.* *.P .real .complex *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#  (From .frc files the appropriatly typed subroutine (REAL/COMPLEX) has to
#  be extracted first.)
#
.f.o:
	$(FC) $(FFLAGS) -c $<
.frc.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >$*_tmp.f
	@echo "$<: $(TYPE) version extracted"
	$(FC) $(FFLAGS) -c $*_tmp.f
	mv $*_tmp.o $*.o
	rm $*_tmp.f
.f.mod:
	$(FC) $(FFLAGS) -c $<
.c.o:
	$(CC) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<
.frc90.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >_tmp_.f90
	@echo "$<: $(TYPE) version extracted"
	$(FC) $(FFLAGS) -c _tmp_.f90
	mv _tmp_.o $*.o
