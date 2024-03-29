include Makefile.in

all : readCharge.so gpoint.so ylm.so yw_excor.so

readCharge.so : readCharge.f90
	$(CMP) $(F2PL) -c $? -m readCharge

gpoint.so : gpoint.f90
	$(CMP) $(F2PL) -c $? -m gpoint

ylm.so : ylm.f90
	$(CMP) $(F2PL) -c $? -m ylm

yw_excor.so : yw_excor.f90
	$(CMP) $(F2PL) -c $? -m yw_excor

clean :
	rm -rf *.so *.o *.a *.mod *.pyc
