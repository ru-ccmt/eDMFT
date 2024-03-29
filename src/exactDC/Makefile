include ../Makefile.in

all : gpoint.so yw_excor.so

gpoint.so : gpoint.f90
	$(CMP) $(F2PL) -c $? -m gpoint
	mv gpoint.*so gpoint.so

yw_excor.so : yw_excor.f90
	$(CMP) $(F2PL) -c $? -m yw_excor
	mv yw_excor.*so yw_excor.so

clean :
	rm -rf *.so *.o *.a *.mod *.pyc
	rm -rf *.dSYM
