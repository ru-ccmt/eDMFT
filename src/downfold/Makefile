include ../Makefile.in

all : latgen.so

latgen.so : latgen.f90
	$(CMP) $(F2PL) -c latgen.f90 -m latgen  --link-lapack_opt
	mv latgen.*so latgen.so

clean :
	rm -f latgen.*so
	rm -rf *.dSYM


