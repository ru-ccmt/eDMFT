# C++ = g++
# CFLAGS = -g
# LIBS = -L/opt/intel/mkl/8.0/lib/32 -lmkl_lapack -lmkl -lguide -lgsl

include ../../Makefile.in

LIBS = $(LLIBS)

CFLAGS = $(OFLAGS) -funroll-all-loops -DNO_ARG_CHECK

all : atomh

atomh : atomh.cc
	$(C++) $(CFLAGS) -o $@ atomh.cc $(LIBS)

clean :
	- rm -f atomh.o atomh
