#PC++ = mpicxx
#C++ = g++
#cc = mpicc
#f77 = ifort
#LLIBS = -L/opt/intel/mkl/8.1/lib/32 -lmkl_lapack -lmkl -lguide
#OFLAGS = -O3 -DNO_ARG_CHECK -funroll-all-loops 
#GFLAGS = -g -DSUNCA_DEBUG
include ../../Makefile.in

LIBS = $(LLIBS)
AOFLAGS = $(OFLAGS) -funroll-all-loops -DNO_ARG_CHECK $(Fmpi_define) # $(OPENMP)
AGFLAGS = $(GFLAGS) -DCTMA_DEBUG $(Fmpi_define)

base_headers = Common.h average.h function.h mesh.h classes.h definitions.h
util_headers = assert.h  complex.h  util.h zeroin.h parser.h  timer.h  blas.h integrate.h

headers = $(base_headers) $(util_headers)

#CFLAGS  $(AGFLAGS)
CFLAGS  = $(AOFLAGS)
FFLAGS = $(AOFLAGS)


all : oca generate

oca : oca.cc
	$(PC++) $(CFLAGS) -o $@ oca.cc $(LIBS)

generate : generate.cc
	$(C++) $(CFLAGS) -o $@ generate.cc


clean :
	- rm -f oca.o oca generate.o generate


.SUFFIXES : .cc
.cc.o:  $(base_headers)
	$(PC++) $(CFLAGS) -c $<

.SUFFIXES : .c
.c.o:  $(base_headers)
	$(pcc) $(CFLAGS) -c $<

.SUFFIXES : .f
.f.o: 
	$(F77) $(FFLAGS) -c $<

