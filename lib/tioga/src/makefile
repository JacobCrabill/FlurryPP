#############################################################################
# Makefile for building TIOGA library file
# Command: make
#############################################################################

####### Compiler, tools and options

MODULENAME=tioga
F90= mpif90
CC = mpicxx
CXX= mpicxx
AR = ar -rvs

# Location of external headers and libraries
MPI_LIB_DIR = -L/usr/local/intel/10.1.011/fce/lib -L/usr/local/openmpi/openmpi-1.4.3/x86_64/ib/intel10/lib

# Location of Flurry headers
FLURRY_INC = ../../../include

# Set Compiler Flags
CFLAGS = -g -std=c++11 -I$(FLURRY_INC) -I.
CXXFLAGS = -g -std=c++11 -I$(FLURRY_INC) -I. -Wno-literal-suffix
FFLAGS = -g -I$(FLURRY_INC)

#INCLUDES = codetypes.h MeshBlock.h ADT.h tioga.h globals.h solver.hpp
#OBJF90 = kaiser.o median.o cellVolume.o
#OBJECTS = buildADTrecursion.o searchADTrecursion.o ADT.o\
#	MeshBlock.o search.o checkContainment.o bookKeeping.o \
#	dataUpdate.o math.o utils.o linklist.o\
#	tioga.o holeMap.o exchangeBoxes.o exchangeSearchData.o exchangeDonors.o\
#	parallelComm.o highOrder.o
#	tiogaInterface.o

OBJF90 = kaiser.o median.o cellVolume.o
OBJECTS = utils.o ADT.o MeshBlock.o parallelComm.o tioga.o

# Set Linker Flags
LDFLAGS= $(MPI_LIB_DIR) -lifcore  -limf -ldl

####### Implicit rules

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CXXFLAGS) -c $< -o $*.o
%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $*.o


####### Build rules

lib:	$(OBJF90) $(OBJECTS) $(INCLUDES)
	$(AR) lib$(MODULENAME).a $(OBJF90) $(OBJECTS)

default: $(OBJF90) $(OBJECTS) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJF90) $(OBJECTS) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -rf *.o lib$(MODULENAME).a

####### Object Compilation

ADT.o: ADT.C ADT.h codetypes.h utils.h
	$(CXX) -c $(CXXFLAGS) -o ADT.o ADT.C

MeshBlock.o: MeshBlock.C MeshBlock.h codetypes.h ADT.h utils.h $(FLURRY_INC)/solver.hpp
	$(CXX) -c $(CXXFLAGS) -o MeshBlock.o MeshBlock.C

parallelComm.o: parallelComm.C parallelComm.h codetypes.h
	$(CXX) -c $(CXXFLAGS) -o parallelComm.o parallelComm.C

utils.o: utils.c codetypes.h utils.h
	$(CXX) -c $(CXXFLAGS) -o utils.o utils.c

tioga.o: tioga.C tioga.h codetypes.h MeshBlock.h parallelComm.h utils.h
	$(CXX) -c $(CXXFLAGS) -o tioga.o tioga.C
		
