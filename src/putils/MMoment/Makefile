CMP = f2py
F2PL = 

all : latgen2.so

latgen2.so : latgen2.f90
	$(CMP) $(F2PL) -c $? -m latgen2
	mv latgen2.*so latgen2.so

clean :
	rm -f latgen2.so *.pyc
	rm -rf *.dSYM
