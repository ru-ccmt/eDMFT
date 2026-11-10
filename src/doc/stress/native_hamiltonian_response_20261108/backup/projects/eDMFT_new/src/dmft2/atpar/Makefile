include ../../Makefile.in

# The discrete radial and complete LAPW/APW+lo tangents/adjoints are opt-in;
# ordinary ATPARN still calls OUTWIN directly. Native response workers can
# reuse these modules without changing the production stress path.
objs = outwin.o  outwin_response.o  lapw_radial_response.o  rint13.o  abc.o  diracout.o  dergl.o  inouh.o  inth.o
#LIBOPT = -dynamiclib

libatpar.a : $(objs) atpar.o readPot.o
	ar  rcs libatpar.a $(objs) atpar.o readPot.o

atpar.so : $(objs) atpar.f90
	$(CMP) -c atpar.f90 -m atpar $(objs)

readPot.so : readPot.f90
	$(CMP) -c readPot.f90 -m readPot

all : libatpar.a atpar.so readPot.so

clean :
	rm -f libatpar.a $(objs) atpar.o atpar.so readPot.so readPot.o *.mod


%.o : %.f90
	$(PF90) $(PIC) $(FFLAGS) -c $< -o $@

# The .mod file must exist before compiling its consumer, also under make -j.
lapw_radial_response.o : outwin_response.o

