include ../../Makefile.in 
#CC = g++-12
#PYBIND = `python3 -m pybind11 --includes` -undefined dynamic_lookup

all : cakw.so

cakw.so : cakw.cc
	$(C++) $(PYBND) -O3 -fopenmp cakw.cc -o cakw.so

clean :
	rm -f cakw.so
