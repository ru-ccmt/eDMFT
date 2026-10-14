include ../../Makefile.in

# C++ = mpic++ # g++
#CFLAGS = -Ofast -std=c++20  -D_MPI # -fsanitize=address -fno-omit-frame-pointer -D_MPI -w
#CFLAGS = -O0 -std=c++20  -D_MPI # -fsanitize=address -fno-omit-frame-pointer -D_MPI -w
#CFLAGS = -g -std=c++20  -fsanitize=address -DBZ_GENERATE_GLOBAL_INSTANCES #-fno-omit-frame-pointer -D_MPI -w 
#CFLAGS = -O2 -fno-fast-math -fno-associative-math -fno-unsafe-math-optimizations -std=c++20 -Wno-deprecated-declarations #-D_MPI # -fsanitize=address -fno-omit-frame-pointer -D_MPI -w
#LIBS = -L/opt/intel/oneapi/mkl/2024.0/lib -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lgsl

#INC = -I/usr/local/include #-I/opt/homebrew/include  -I/usr/include 
#LIBS = -framework Accelerate -L/opt/homebrew/lib -lgsl

CFLAGS = $(PFLAGS) $(Fmpi_define) $(GSLINC) -I. -std=c++20
LIBS = $(PLIBS)

all : bctqmc

bctqmc : boldc.cc readinput.h interpolate.h measureW.h skeletons.h wfile.h basics.h generatediags.h readinput.h tanmesh.h util.h wdiags.h wevaluate.h whybrid.h wsign.h Gexp.h dynamic_mesh.h any_axis.h real_axis.h imag_axis.h nca_oca_imag_axis.h nca_oca_real_axis.h imag_axis_spectra.h real_axis_spectra.h frequency_loops_real_axis.h wbroad.h read_cix.h wmpi.h
	$(PC++) $(CFLAGS) boldc.cc -o $@ $(LIBS)

clean:
	rm -f bctqmc boldc.o
