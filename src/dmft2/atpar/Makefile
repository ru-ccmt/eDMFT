include ../../Makefile.in

# The discrete radial and complete LAPW/APW+lo tangents/adjoints are opt-in;
# ordinary ATPARN still calls OUTWIN directly. Native response workers can
# reuse these modules without changing the production stress path.
objs = outwin.o  outwin_response.o  lapw_radial_response.o  lapw_radial_hamiltonian.o  lapw_radial_potential.o  rint13.o  abc.o  diracout.o  dergl.o  inouh.o  inth.o
#LIBOPT = -dynamiclib

libatpar.a : $(objs) atpar.o readPot.o
	ar  rcs libatpar.a $(objs) atpar.o readPot.o

atpar.so : $(objs) atpar.f90
	$(CMP) -c atpar.f90 -m atpar $(objs)

readPot.so : readPot.f90
	$(CMP) -c readPot.f90 -m readPot

all : libatpar.a atpar.so readPot.so

clean :
	rm -f libatpar.a $(objs) atpar.o atpar.so readPot.so readPot.o *.mod lapw_radial_capi.o liblapw_response.so


%.o : %.f90
	$(PF90) $(PIC) $(FFLAGS) -c $< -o $@

# The .mod file must exist before compiling its consumer, also under make -j.
lapw_radial_response.o : outwin_response.o

lapw_radial_hamiltonian.o : lapw_radial_response.o

# Mixed-l nonspherical integrals share the primitive radial representation.
lapw_radial_potential.o : lapw_radial_response.o

# Optional in-memory response API. This is not part of the ordinary ATPARN
# archive or the default dmft2 path. Prepare radial contexts once, reuse them
# across k points and in the forward/adjoint solves. Override the linker flag
# if a platform requires a different shared-library option.
RADIAL_SHARED_LDFLAGS ?= -shared
.PHONY: response
response : liblapw_response.so

liblapw_response.so : outwin.o outwin_response.o lapw_radial_response.o lapw_radial_capi.o
	$(PF90) $(PIC) $(FFLAGS) $(RADIAL_SHARED_LDFLAGS) -o $@ $^

lapw_radial_capi.o : lapw_radial_response.o
