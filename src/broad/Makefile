include ../Makefile.in

# C++ = g++
# CFLAGS = -Ofast -std=c++2a

LIBS = $(LALIB) $(GSLLIB) -I../includes # -framework Accelerate -lgsl

all : broad

#broad : broad.cc
#	$(C++) $(PFLAGS) -o $@ $? $(LIBS)

broad : nbrd.cc
	$(C++) $(PFLAGS) $(OFLAGS) -o $@ $? $(LIBS)

clean :
	rm -f broad

