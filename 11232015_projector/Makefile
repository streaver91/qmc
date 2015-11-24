FC=gfortran
FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -fbounds-check -fbacktrace
F90FLAGS = -finit-local-zero -O0 -ffree-line-length-none -x f95-cpp-input -Wall -g -fbounds-check -fbacktrace
#FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -Wall
#F90FLAGS = -finit-local-zero -O0 -ffree-line-length-none -x f95-cpp-input -Wall
QUENCH  = -L$(HOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
LINPACK = -L$(HOME)/lib2/linpack -llinpack
BLAS = -L$(HOME)/lib2/blas -lblas
#LIBS = $(LINPACK) $(QUENCH)
#LIBS = $(QUENCH) $(LINPACK) $(BLAS)
LIBS = $(BLAS)
LIBS = -llapack-3
#LIBS = -L/usr/lib/libblas -lblas -L/usr/lib -llapack_atlas
LIBS = -L/home/cyrus/lib2/lapack -llapack -L/usr/lib/libblas -lblas -L/home/cyrus/lib -lcyrus

.SUFFIXES:
.SUFFIXES: .f90 .f95 .o .f .c

.f90.o:
	$(FC) $(F90FLAGS) -o $@ -c $<

OBJS = types.o rannyu.o tools.o more_tools.o fixed_node.o matinv.o

all: fixed_node
#all: fixed_node fixed_node2

fixed_node: $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@
