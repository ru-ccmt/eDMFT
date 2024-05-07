include ../../Makefile.in

BLITZ = -I../../includes

all : gaunt.so dpybind.so

gaunt.so: gaunt.f90
	$(CMP) $(F2PL) -c $? -m gaunt 
	mv gaunt.*so gaunt.so

dpybind.so  : dpybind.cc
	$(C++) $(BLITZ) $(PYBND) -O2 dpybind.cc -o dpybind.so

clean:
	rm -rf build gaunt.*so dpybind.*.so dpybind.so
	rm -rf *.dSYM 

