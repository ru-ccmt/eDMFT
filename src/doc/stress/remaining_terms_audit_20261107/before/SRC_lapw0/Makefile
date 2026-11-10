###############################################################################
#                                                                             #
#                          Generic Makefile for lapw0                         #
#                                                                             #
#  make           ... generate executable for the sequential version          #
#  make seq       ... generate executable for the sequential version          #
#  make para      ... generate executable for the parallel version            #
#  make all       ... generate sequential and parallel version                #
#  make last      ... repeat latest make (seq or para)                        #
#  make clean     ... delete unnecessary files                                #
#                                                                             #
###############################################################################

.SUFFIXES:        .F
.SUFFIXES:        .F90
SHELL = /bin/sh

###############################################################################
#                           User defined variables                            #
###############################################################################

#---------------- compilers ----------------
FC = gfortran
MPF = 
CC = gcc

#------------------ libxc ------------------
LIBXCROOT = 
LIBXC_FORTRAN = 
LIBXC_LIBDNAME = 
LIBXC_LIBNAME = 
ifneq ($(LIBXCROOT),)
  LIBXC_FOPT = -DLIBXC -I$(LIBXCROOT)include 
  LIBXC_LIBS = -L$(LIBXCROOT)$(LIBXC_LIBDNAME) -l$(LIBXC_FORTRAN) -l$(LIBXC_LIBNAME)
endif

#------------------- fftw ------------------
FFTW_VERSION = FFTW3
FFTWROOT = /opt/homebrew/
FFTW_LIB = lib
FFTW_LIBNAME = fftw3
ifneq ($(FFTWROOT),)
  FFTW_FOPT = -D$(FFTW_VERSION) -DFFTW_OMP -I$(FFTWROOT)include 
  FFTW_LIBS = -L$(FFTWROOT)$(FFTW_LIB) -l$(FFTW_LIBNAME)_omp -l$(FFTW_LIBNAME)
  FFTW_PLIBS = -l$(FFTW_LIBNAME)_mpi 
endif

#------------- compiler options ------------
OMP_SWITCH = -fopenmp
FOPT =  -ffree-form -O2 -ftree-vectorize -march=native -ffree-line-length-none -ffpe-summary=none $(LIBXC_FOPT) $(FFTW_FOPT) $(OMP_SWITCH)
# Diagnostic alternative: replace -O2 with -Og -g -fcheck=all -fbacktrace.
FPOPT =  -O -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I$(MKLROOT)/include $(LIBXC_FOPT) $(FFTW_FOPT) $(OMP_SWITCH)
DParallel = 
FGEN = $(PARALLEL)

#-------------- loader options -------------
LDFLAGS = $(FOPT) -L../SRC_lib
R_LIBS     = $(FFTW_LIBS) -framework Accelerate $(LIBXC_LIBS) 
RP_LIBS = $(FFTW_PLIBS) $(R_LIBS)

###############################################################################
#                       End of user defined variables                         #
###############################################################################

DESTDIR  = ./

S_EXECNAME = lapw0
P_EXECNAME = lapw0_mpi

###############################################################################
###############################################################################

S_EXEC   = $(DESTDIR)$(S_EXECNAME)
P_EXEC   = $(DESTDIR)$(P_EXECNAME)

#..............................................................................
#
#  All object files, including $(SP_OBJS) and $(OBJSMOD) 
#
OBJS  =cputim.o modules.o libxc_mod.o W2kinit.o fft_modules.o reallocate.o \
        acggac.o ainv.o ak13.o am05_xscss.o b88.o br89.o brj.o cap.o charge.o \
        chfac.o chslv.o corgga.o corpbe_revtpss.o corpbe_tpss.o \
        cub_xc_back.o corlsd.o dergl.o dfxhpbe.o dfxrevtpss.o dfxtpss.o dfundr.o \
        dylm.o efg.o energy.o epot1.o eramps.o errclr.o errflg.o ev92.o \
        ev92ex.o exch.o fxhpbe.o \
        fx_revtpss.o fx_tpss.o gbass.o gbj.o gcor.o getff1.o getfft.o \
        gpoint.o gpointL.o gpointm.o grans.o gtfnam.o hcthc.o hcthx.o \
        hjs.o htbs.o id_func_get.o ifflim.o inputpars.o kcis.o \
        lapw0.o latgen.o lb94.o libxc.o lyp.o mbeefx.o mgga_ms.o multfc.o multsu.o mvsx.o \
        optb86b.o optb88.o outerr.o pbea.o pbec.o pbeint.o pbex.o poissn.o \
        potfac.o pwxad4.o pwxad5.o readstruct.o rean_direct_fft.o rean0.o rean3.o \
        rean4.o rge2.o rhopw.o rotate.o rotdef.o rpbe.o rscanc.o rscanx.o r2scanc.o r2scanx.o \
        scanc.o scanx.o setff0.o setff1.o setfft.o setff2.o sortag.o \
        sg4c.o sg4x.o sicpbe.o sicpbe_revtpss.o sicpbe_tm.o sicpbe_tpss.o sphbes.o \
        spline.o srolyl.o stern.o sumfac.o suml.o SymmRot.o tmx.o vpw91.o \
        vresp.o vs98.o vwn5.o vxclm2.o \
        vxcpw2.o vx_screened.o wc05.o workf1.o W2kutils.o xcener.o \
        xcpot.o xcpot1.o xcpot1Q.o xcpot3.o ykav.o ylm.o coulint.o \
        c_alpha_m.o gaunt.o notri.o t3j0.o t3j.o ph.o finl_elect_str.o \
        kahan_summ_gxyz.o stress_gga.o
#  Object file for parallel and sequential
SP_OBJS = modules.o libxc_mod.o W2kinit.o fft_modules.o reallocate.o energy.o \
          finl_elect_str.o getff1.o getfft.o gtfnam.o lapw0.o outerr.o \
          rean_direct_fft.o rean0.o rean3.o rean4.o setff1.o setff2.o setfft.o \
          stress_gga.o xcpot1.o xcpot1Q.o xcpot3.o eramps.o inputpars.o

#  Object files using modules
OBJSMOD = br89.o brj.o charge.o dfundr.o energy.o epot1.o eramps.o finl_elect_str.o \
          gbj.o getff1.o gpoint.o gpointL.o gtfnam.o hcthc.o hcthx.o hjs.o \
          inputpars.o lapw0.o latgen.o libxc.o mgga_ms.o multsu.o outerr.o pbea.o pbec.o \
          pbex.o poissn.o pwxad4.o pwxad5.o readstruct.o rean_direct_fft.o rean0.o \
          rean3.o rean4.o rotdef.o rpbe.o setff1.o setff2.o setfft.o stern.o stress_gga.o \
          vresp.o vxclm2.o vx_screened.o xcpot1.o xcpot1Q.o xcpot3.o ykav.o  
#........................................................................
default: seq

all:    
	make seq para

last: 
	if [ ! -f .sequential -a ! -f .parallel ]; then make seq; fi
	if [ -f .sequential ]; then make seq; fi
	if [ -f .parallel ]; then make para; fi
        
seq: keep_s_files
	$(MAKE) $(S_EXEC) FORT=$(FC) FFLAGS='$(FGEN) $(FOPT)'
	if [ -f $(LIBXCROOT)/include/xc_funcs.h ]; then cp $(LIBXCROOT)/include/xc_funcs.h .; fi

para: keep_p_files
	$(MAKE) PARALLEL=$(DParallel) $(P_EXEC) \
	  FORT=$(MPF) FFLAGS='$(FGEN) $(FPOPT) $(DParallel)'
	if [ -f $(LIBXCROOT)/include/xc_funcs.h ]; then cp $(LIBXCROOT)/include/xc_funcs.h .; fi

#..............................................................................

$(S_EXEC): $(OBJS)
	$(FORT) -o $(S_EXEC) $(OBJS) $(LDFLAGS) $(R_LIBS)

$(P_EXEC): $(OBJS)
	$(FORT) -o $(P_EXEC) $(OBJS) $(LDFLAGS)  $(RP_LIBS) 

# sequential
keep_s_files:
	if [ -f .parallel ]; then \
	   rm -f .parallel $(SP_OBJS) *.mod; \
	fi
	touch .sequential

# parallel
keep_p_files:
	if [ -f .sequential ]; then \
	   rm -f .sequential $(SP_OBJS) *.mod; \
	fi
	touch .parallel

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):		param.inc Makefile
$(OBJSMOD):		modules.o fft_modules.o W2kinit.o pwxad4.o


#..............................................................................
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#
.f.o:
	$(FORT) $(FFLAGS) -c $<
.F.o:
	$(FORT) $(FFLAGS) -c $<
.c.o:
	$(CC) -c $<

.f90.o:
	$(FORT) $(FFLAGS) -c $<
.F90.o:
	$(FORT) $(FFLAGS) -c $<
