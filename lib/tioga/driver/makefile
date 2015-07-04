MODULENAME=tioga
F90= mpif90
CC = mpicc
CXX= mpicxx
AR = ar -rvs
FFLAGS= -O3 #-g -fdefault-real-8 -fcheck=bounds
#FFLAGS=-CB -traceback -r8
CFLAGS = -O3 #-g -traceback #-fpe0
INCLUDES = 
OBJEXEC=
OBJECTS=
OBJF90 = modules.o gridutils.o testTioga.o
TIOGALIB=../src.d/libtioga.a

#LDFLAGS= -L../src/ -L/usr/local/intel/10.1.011/fce/lib -L/usr/local/openmpi/openmpi-1.4.3/x86_64/ib/intel10/lib  -lifcore  -limf -ldl  -ltioga -lstdc++

LDFLAGS= -L../src/ -ltioga -lstdc++ -lmpi_cxx

default: $(OBJF90) $(LIB)
	$(F90) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -rf *.o $(MODULENAME).exe

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
