# C++ = g++
# CFLAGS = -g
# LIBS = 
include ../../Makefile.in
LIBS = 
CFLAGS = $(OFLAGS) -funroll-all-loops -DNO_ARG_CHECK

all : skrams

skrams : skrams.cc
	$(C++) $(CFLAGS) -o skrams skrams.cc $(LIBS)


clean :
	- rm -f skrams.o skrams
