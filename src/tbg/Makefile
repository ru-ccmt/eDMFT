include ../Makefile.in

FLINK = $(F2PY_LAPACK) -lgomp
FLAGS = --opt='-O3 -fopenmp'

#FLAGS = -L/opt/intel/lib/intel64/ -liomp5 -lpthread --opt='-O2 -fopenmp'  --link-lapack_opt --fcompiler=intelem
#C++ = g++-14
#BLITZ = -I. 
#PYBND = `python -m pybind11 --includes` -shared -std=c++11 -fPIC -undefined dynamic_lookup
#FLINK = -lgomp --link-lapack_opt


all : tbg_scc.so

tbg_scc.so : tbg_scc.f90
	f2py $(FLAGS) -c tbg_scc.f90 -m tbg_scc $(FLINK)
	mv tbg_scc.*so tbg_scc.so

clean : 
	rm -f tbg_scc.so *.cpython-*.so
