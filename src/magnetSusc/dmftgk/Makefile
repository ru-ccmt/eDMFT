include ../../Makefile.in

EXEC = dmftgk
OBJS = dmftgk.o eigsys.o selfenergy.o zinv.o

all : $(EXEC)

$(EXEC): $(OBJS)
	$(PF90) $(FFLAGS) -o $(EXEC) $(OBJS) $(WLIBS)


clean :
	rm -f $(OBJS) $(EXEC) *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#
.f.o:
	$(PF90) $(FFLAGS) $(debug) -c $<
.c.o:
	$(CC)  -c $<


.SUFFIXES : .f90
.f90.o:
	$(PF90) $(FFLAGS) $(debug) -c $<

.SUFFIXES : .F
.F.o:
	$(preproc) -P  -C -traditional $(Fmpi_define) $< $*.f90
	$(PF90) $(FFLAGS) $(debug) -c $*.f90
	rm $*.f90


