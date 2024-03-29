include ../../Makefile.in
#INSTALL = example_FeO

all : chi2f.so

chi2f.so : chi2f.f90
	$(CMP) -m chi2f -c chi2f.f90
	mv chi2f.*so chi2f.so
clean : 
	rm -f chi2f.*so
	rm -rf *.dSYM*
install:
	cp pmesh.py ancont.py chi2f.*so strans run.sh $(INSTALL)
