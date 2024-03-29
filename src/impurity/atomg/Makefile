# C++ = g++
# OFLAGS = -O3 -DNO_ARG_CHECK
# GFLAGS = -g -DCTMA_DEBUG 
# LIBS = -L/opt/intel/mkl/8.1/lib/32 -lmkl_lapack -lmkl -lguide

include ../../Makefile.in
AOFLAGS = $(OFLAGS)
AGFLAGS = $(GFLAGS)
LIBS = $(LLIBS)

headers = assert.h  blas.h complex.h  function.h mesh.h util.h parser.h

#CFLAGS = $(AGFLAGS)
CFLAGS = $(AOFLAGS)

all : atom

atom : atom.cc $(headers)
	$(C++) $(CFLAGS) -o $@ atom.cc $(LIBS)

clean :
	- rm -f atom
