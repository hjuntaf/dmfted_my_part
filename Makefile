CC = mpic++

RM = rm -f 

LDFLAGS =  -fopenmp -mcmodel=medium -L/opt/lapack/3.8.0/gcc-4.8.5/lib -L/opt/lapack/3.8.0/gcc-4.8.5/lib64  -llapack -L/opt/atlas/3.10.3/gcc-4.8.5/lib/  -latlas -L/opt/gsl/2.5/gcc-4.8.5/lib -L/opt/mpi/gcc-4.8.5/openmpi-2.1.5/lib64 -L/home/arago/local/arpackpp-source/external  -lgfortran -lgsl -lgslcblas  -larpack -lsuperlu 
#LDFLAGS = -fopenmp -mcmodel=medium -L/opt/lapack/3.8.0/gcc-4.8.5/lib/ -L/opt/blas/3.8.0/gcc-4.8.5/lib/ -L/opt/gsl/2.5/gcc-4.8.5/lib -L/opt/mpi/gcc-4.8.5/openmpi-2.1.5/lib64 -L/home/arago/local/arpackpp-source/external  -lgfortran -lgsl -lgslcblas -llapack -lblas -larpack -lsuperlu
#LDFLAGS = -fopenmp -mcmodel=medium -L/scratch/usr/local/gsl2.1/lib -lgsl -lgslcblas -lblas -larpack -llapack -lsuperlu -L/opt/openmpi-1.10.3/lib64 -L/home/arago/local/arpackpp-source/external 

#CFLAGS = -fgnu89-inline -c -g -O2 -Wall -mcmodel=medium -fopenmp -I${HOME}/local/include/arpackpp
#CFLAGS = -c -g -O2 -std=c99 -Wall -mcmodel=medium -fopenmp -I${HOME}/local/include/arpackpp
CFLAGS = -c -g -O2 -Wall -Wno-restrict -Wno-maybe-uninitialized -Wno-misleading-indentation -Wno-format-overflow -Wno-write-strings -mcmodel=medium -fopenmp  -I${HOME}/local/include/arpackpp -I/opt/gsl/2.5/gcc-4.8.5/include/ 
BUILD = run vert local static
OBJS = run.o matrix.o basic.o hamil.o ls_basis.o arp.o

.PHONY: all clean  

default: run
all: $(BUILD) 
clean: ; $(RM) *.o 

run: $(OBJS)
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o $@
vert: vert_new.o matrix.o basic.o
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o vert_new
local: local.o matrix.o basic.o ls_basis.o hamil.o
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o $@
localwsq: local_wsqbroad.o matrix.o basic.o ls_basis.o hamil.o
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o $@

static: static.o matrix.o basic.o operators.o
	$(CC) $^    $(LDFLAGS) $(LDLIBS) -o $@
moment: moment.o matrix.o basic.o operators.o
	$(CC) $^    $(LDFLAGS) $(LDLIBS) -o $@


chi: chi.o matrix.o basic.o hamil.o ls_basis.o operators.o
	$(CC) $^    $(LDFLAGS) $(LDLIBS) -o $@

run.o:	$($@:.o=.c)	lattice.h matrix.h matrix.c inverse.c ls_basis.h ls_basis.c basic.c ls.c hop.c memcifbathm.c ran2.c hamil.h hamil.c arp.h arp.c Makefile 
chi.o:		$($@:.o=.c)	lattice.h matrix.h matrix.c inverse.c ls_basis.h ls_basis.c basic.c ls.c hop.c memcifbathm.c ran2.c hamil.h hamil.c arp.h arp.c response.c Makefile 
ls_basis.o:	$($@:.o=.c)	lattice.h matrix.h matrix.c ls_basis.h ls_basis.c Makefile 
arp.o:		$($@:.o=.c)	Makefile
hamil.o:	$($@:.o=.c)	Makefile
basic.o:	$($@:.o=.c)	Makefile
operators.o:	$($@:.o=.c)	Makefile
matrix.o:	$($@:.o=.c)	matrix.h	Makefile
static.o:	$($@:.o=.c)	ls_basis.c basic.h basic.c lattice.h operators.c	Makefile
moment.o:	$($@:.o=.c)	ls_basis.c basic.h basic.c lattice.h operators.c	Makefile
vert_new.o:	$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile
local.o:	$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile
local_wsqbroad.o:	$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile

.c.o:
	$(CC) $(CFLAGS) -o $@ $< 
