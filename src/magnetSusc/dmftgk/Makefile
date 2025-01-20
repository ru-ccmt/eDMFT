include ../../Makefile.in

EXEC = dmftgk
OBJS = mod_mpi.o dmftgk.o eigsys.o selfenergy.o zinv.o

all : $(EXEC)

$(EXEC): $(OBJS)
	$(PF90) $(WFOPT) -o $(EXEC) $(OBJS) $(WLDFLAGS) $(WLIBS) 

clean :
	rm -f $(OBJS) $(EXEC) *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#
%.o : %.f90
	$(PF90) $(WFOPT) -c $< -o $@

%.o : %.F90
	$(PF90) $(WFOPT) -c $< -o $@

%.o : %.FPP
	$(PF90) $(Fmpi_define) $(WFOPT) -c $< -o $@

%.o : %.c
	$(CC) -c $< -o $@

%.so : %.f90
	$(CMP) -c $< -m $*
