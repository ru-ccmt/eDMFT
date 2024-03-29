include ../../Makefile.in

all : w2k_nn.so

w2k_nn.so : w2k_nn.f90
	$(CMP) $(F2PL) -c $? -m w2k_nn
	mv w2k_nn.*so w2k_nn.so

clean :
	rm -rf *.so *.o *.a *.mod
	rm -rf *.dSYM
