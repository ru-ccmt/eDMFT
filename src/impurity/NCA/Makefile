# C++ = mpicxx
# cc = mpicc
# f77 = ifort
# LIBS = -L/opt/intel/mkl/8.1/lib/32 -lmkl_lapack -lmkl -lguide
# OFLAGS = -O3 -DNO_ARG_CHECK -funroll-all-loops
# GFLAGS = -g -DSUNCA_DEBUG 

include ../../Makefile.in

AOFLAGS = $(OFLAGS) -funroll-all-loops -DNO_ARG_CHECK $(OPENMP)
AGFLAGS = $(GFLAGS) -DCTMA_DEBUG
LIBS = $(LLIBS)

base_headers = Common.h average.h function.h mesh.h
util_headers = assert.h  complex.h  util.h zeroin.h parser.h timer.h

headers = $(base_headers) $(util_headers)

CFLAGS  = $(AOFLAGS) # $(AGFLAGS)

all : nca generate

nca : nca.cc Common.h
	$(C++) $(CFLAGS) -o $@ nca.cc $(LIBS)

generate : generate.cc
	$(C++) $(CFLAGS) -o generate generate.cc

clean :
	- rm -f nca.o generate.o nca generate 


.SUFFIXES : .cc
.cc.o:  $(base_headers)
	$(C++) $(CFLAGS) -c $<

.SUFFIXES : .c
.c.o:  $(base_headers)
	$(cc) $(cFLAGS) -c $<

.SUFFIXES : .f
.f.o: 
	$(F77) $(FFLAGS) -c $<

