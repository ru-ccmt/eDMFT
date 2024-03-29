include ../Makefile.in

FORT = $(PF90)
LFLIBS = $(WLDFLAGS) $(WLIBS) -L./lattice -L./atpar -llattice -latpar
FFLAGS = $(WFOPT)  -I./lattice -I./atpar
LDFLAGS = $(WFOPT)

OBJS = mpi.o modules.o abclm.o couple.o errclr.o  errflg.o garadme.o gtfnam.o  \
	harmon2.o  atpar.o    read_inso.o  find_nlos.o   \
        lapwso.o   cali.o    sphbes.o vnsrint.o  \
        lomain.o        cputim.o  dvbes1.o     lagder.o  \
        ylm.o      sphbru.o   outerr.o  lap_bp.o \
        inouh.o    diracout.o inth.o    select.o   abc_r.o   rlomain.o \
        kptout.o angle.o gaunt1.o gaunt2.o getmaxdim.o \
        dergl.o    hscalc.o   dergl2.o  vderiv.o  read_def_file.o read_orbital.o \
        symop.o garadorb.o hns.o hsocalc.o horbcalc.o vorblo.o euler.o hmsec.o get_nloat.o get_new_nloat.o

all : llattice latpar lapwso 

llattice:
	cd lattice && ${MAKE} liblattice.a

latpar :
	cd atpar && ${MAKE} libatpar.a

lapwso : $(OBJS) llattice latpar
	$(FORT) -o $@ $(OBJS) $(LDFLAGS) $(LFLIBS)


clean:
	-rm -f *.o _tmp_.* *.P .real .complex *.mod *__genmod.f90 lapwso
	cd lattice; make clean; cd ../
	cd atpar; make clean; cd ../

#...................
# Pattern rules
#...................
%.o : %.f90
	$(FORT) $(FFLAGS) -c $< -o $@

%.o : %.FPP
	$(FORT) $(FFLAGS) $(Fmpi_define) -c $< -o $@

%.o : %.c
	$(CC) -c $< -o $@
