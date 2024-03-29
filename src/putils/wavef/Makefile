include ../../Makefile.in

all : w2k_atpar.so

w2k_atpar.so : w2k_atpar.f90
	$(CMP) $(F2PL) -c $? -m w2k_atpar
	mv w2k_atpar.*so w2k_atpar.so

clean :
	rm -rf *.so *.o *.a *.mod
	rm -rf *.dSYM
