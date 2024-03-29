include ../../Makefile.in

objs = structure.o inverssymdef.o generate_lattice.o correct_rotij.o get_rotij_tauij.o

liblattice.a : $(objs)
	ar rcs liblattice.a $(objs)

clean :
	rm -f $(objs) liblattice.a *.mod

%.o : %.f90
	$(PF90) $(PIC) $(FFLAGS) -c $< -o $@

