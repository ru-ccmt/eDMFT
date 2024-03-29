include ../../Makefile.in

all : libclebsch.a krams.so strns.so linalg.so ldau.so brd.so

libinalg.a : libinalg.f90 lib_cinv.f
	$(F90) -c $(PIC) $(FFLAGS) libinalg.f90 
	$(F77) -c $(PIC) $(FFLAGS) lib_cinv.f
	ar cr libinalg.a libinalg.o lib_cinv.o

libminim.a : minim.f
	$(F77) -c $(PIC) $(FFLAGS) minim.f
	ar cr libminim.a minim.o

libSk.a : mod_Sk.f90
	$(F90) -c $(PIC) $(FFLAGS) $?
	ar cr $@ mod_Sk.o

libclebsch.a : mod_clebsch.f90
	$(F90) -c $(PIC) $(FFLAGS) mod_clebsch.f90
	ar cr libclebsch.a mod_clebsch.o

fort_read.so : fort_read.f90
	$(CMP) $(F2PL) -c $? -m fort_read
	mv fort_read.*so fort_read.so

mod_Sk.so : mod_Sk.f90
	$(CMP) $(F2PL) -c $? -m mod_Sk
	mv mod_Sk.*so mod_Sk.so

symmetr.so : symmetr.f90 rot_mat0.f90 rot_mat1.f90
	$(CMP) $(F2PL) -c symmetr.f90 rot_mat1.f90 rot_mat0.f90 -m symmetr $(F2PY_LAPACK) -L. -linalg
	mv symmetr.*so symmetr.so

strns.so : strns.f90
	$(CMP) $(F2PL) -c $? -m strns -L. -lclebsch
	mv strns.*so strns.so

transb.so : transb.f90
	$(CMP) $(F2PL) -c $? -m transb  $(F2PY_LAPACK) -L. -lSk -linalg -lminim
	mv transb.*so transb.so

cmp_olocef.so : cmp_olocef.f90
	$(CMP) $(F2PL) -c $? -m cmp_olocef $(F2PY_LAPACK) -L. -linalg
	mv cmp_olocef.*so cmp_olocef.so

linalg.so : libinalg.f90
	$(CMP) $(F2PL) -c $? -m linalg $(F2PY_LAPACK)
	mv linalg.*so linalg.so

tetra.so : wtetrah.f90 tetrah.f90 libtetrah.f90 libtetrahi.f90 creorder.f90 libmetrop.f90 tetrahi.f90 sort.f
	$(F77) -c $(FFLAGS) sort.f
	$(F90) $(FFLAGS) -c libtetrah.f90
	$(F90) $(FFLAGS) -c libtetrahi.f90
	$(F90) $(FFLAGS) -c creorder.f90
	ar cr libtetrah.a libtetrah.o libtetrahi.o creorder.o sort.o
	$(F90) $(FFLAGS) -c libmetrop.f90
	ar cr libmetrop.a libmetrop.o
	$(CMP) $(F2PL) -c wtetrah.f90 tetrah.f90 tetrahi.f90 -m tetra -L. -ltetrah -lmetrop
	mv tetra.*so tetra.so

cmp_rho.so : cmp_rho1.f90
	$(CMP) $(F2PL) -c cmp_rho1.f90 -m cmp_rho $(F2PY_LAPACK) -L. -lSk -linalg
	mv cmp_rho.*so cmp_rho.so

mix.so : mix.f90
	$(CMP) $(F2PL) -c mix.f90 -m mix $(F2PY_LAPACK)
	mv mix.*so mix.so 

ldau.so : ldau.f90
	$(CMP) $(F2PL) -c ldau.f90 -m ldau -L. -lclebsch
	mv ldau.*so ldau.so

brd.so : brd.f90
	$(CMP) $(F2PL) -c brd.f90 -m brd
	mv brd.*so brd.so 

dia_gho.so : dia_gho.f90
	$(CMP) $(F2PL) -c dia_gho.f90 -m dia_gho  $(F2PY_LAPACK)
	mv dia_gho.*so dia_gho.so

mtoverlap.so : mtoverlap.f90
	$(CMP) $(F2PL) -c $? -m mtoverlap  $(F2PY_LAPACK) -L. -lSk -linalg
	mv mtoverlap.*so mtoverlap.so

krams.so : krams.f90
	$(CMP) $(F2PL) -c $? -m krams
	mv krams.*so krams.so

smf.so : smf.f90
	$(CMP) $(F2PL) -c smf.f90 -m smf $(F2PY_LAPACK)
	mv smf.*so smf.so

optics.so : optics.f90
	$(CMP) $(F2PL) -c optics.f90 dia_gho.f90 -m optics  -L. -linalg $(F2PY_LAPACK)
	mv optics.*so optics.so

libintgs.a : intgs.f90 Slatec.f
	$(F77) -c $(PIC) $(FFLAGS) Slatec.f
	$(F90) -c $(PIC) $(FFLAGS) intgs.f90
	ar cr libintgs.a intgs.o Slatec.o

cmtrans.so : cmtrans.f90 intgs.f90 libintgs.a
	$(CMP) $(F2PL) -c cmtrans.f90 dia_gho.f90 -m cmtrans -L. -lintgs -linalg $(F2PY_LAPACK)
	mv cmtrans.*so cmtrans.so

YlmAYlm.so : YlmAYlm.f90
	$(CMP) $(F2PL) -c YlmAYlm.f90 -m YlmAYlm
	mv YlmAYlm.*so YlmAYlm.so

libYlm.so : libYlm.f
	$(CMP) $(F2PL) -c libYlm.f -m libYlm
	mv libYlm.*so libYlm.so

clean :
	rm -rf *.so *.o *.a *.mod
	rm -rf *.dSYM
