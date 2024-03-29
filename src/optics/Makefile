include ../Makefile.in

EXEC = dmftopt
OBJS = mod_mpi.o optmain.o eigsys.o selfenergy.o interpolate.o


all : $(EXEC)
debug = #-check all

$(EXEC): $(OBJS)
	$(PF90) $(FFLAGS) $(debug) -o $(EXEC) $(OBJS) $(WLDFLAGS) $(WLIBS)


clean :
	rm -f $(OBJS) $(EXEC) *.mod mod_mpi.f90

#...................
# Pattern rules
#...................
%.o : %.f90
	$(PF90) $(FFLAGS) $(debug) -c $< -o $@

%.o : %.c
	$(CC) -c $< -o $@

%.o : %.FPP
	$(PF90) $(Fmpi_define) $(FFLAGS) $(debug) -c $< -o $@

