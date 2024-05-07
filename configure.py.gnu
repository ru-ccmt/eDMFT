#! /usr/bin/env python
# -*- coding: utf-8 -*-

class Config:
  prefix      = "bin"    # Installation path

  compiler    = "GNU"          # Compiler
  fc          = "gfortran"     # Fortran compiler
  cc          = "gcc"          # C compiler
  cxx         = "g++"          # C++ compiler


  cflags      = "-O2"          # linker flags for C programs
  fflags      = "-O2 -fallow-argument-mismatch"          # linker flags for Fortran programs
  ldflags      = ""             # linker flags debuggin programs
  ompflag     = "-fopenmp"     # linker/compiler flag for openmp

  mpi_define  = "-D_MPI"       #
  pcc         = "mpicc"        # C compiler 
  pcxx        = "mpicxx"       # C++ compiler 
  pfc         = "mpif90"       # Fortran compiler 
  
  blasname    = "GNU"          # BLAS   library
  blaslib     = "-L/opt/gnu/12.3.0/OpenBLAS/lib -lopenblas"         # BLAS   library
  lapacklib   = "-L/opt/gnu/12.3.0/lib64 -llapack"             # LAPACK library
  fftwlib     = "-lfftw3_omp -lfftw3"  # FFTW   library
  gsl         = "-lgslcblas -lgsl"     # GSL    library


  f2pylib     = "--f90flags='-openmp '"	       # F2PY   library	
  f2pyflag    = "--opt='-O2'"	       # F2PY   library	

  ranlib      = ""             # Ranlib
  arflags     = "rc"           # ar flags

  make        = "make"
  def __init__(self, version):
    self.version = version
  def __getattr__(self,key):
    return None
