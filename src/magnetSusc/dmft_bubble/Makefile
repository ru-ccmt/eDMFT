include ../../Makefile.in

#CFFLAGS = $(FFLAGS) -check all

OBJS = utils.o linlogm.o klist.o greenk.o cputim.o dmft_real_bubble.o

all : dmft_real_bubble

dmft_real_bubble : $(OBJS)
	$(PF90) $(FFLAGS) -o dmft_real_bubble $(OBJS) $(WLIBS)


clean:
	rm -f $(OBJS) dmft_real_bubble *.mod

.SUFFIXES : .c
.c.o:
	$(CC)  -c $<

.SUFFIXES : .f90
.f90.o:
	$(PF90) $(FFLAGS) -c $<
