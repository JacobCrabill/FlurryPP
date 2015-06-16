#############################################################################
# Makefile for building Flurry
# Command: make debug
#          make release
#          make openmp
#############################################################################

####### Compiler, tools and options

CXX           = g++
LINK          = g++
MPICXX        = mpicxx
MPILD         = mpicxx
LIBS          = $(SUBLIBS)

# Location of libmetis.a, metis.h
METIS_LIB_DIR = /usr/lib/
METIS_INC_DIR = /usr/include/

# Location of mpi.h 
MPI_INC_DIR   = /usr/lib/openmpi/include

CXX_BASE    = -pipe -Wall -W -std=c++11 $(DEFINES)
CXX_STD     = -g -02
CXX_DEBUG   = -g -pg -O0
CXX_RELEASE = -O3

CXXFLAGS_RELEASE = $(CXX_BASE) $(CXX_RELEASE) -D_NO_MPI $(DEFINES)
CXXFLAGS_DEBUG   = $(CXX_BASE) $(CXX_DEBUG) -D_NO_MPI $(DEFINES)
CXXFLAGS_OPENMP  = $(CXX_BASE) $(CXX_RELEASE) -fopenmp -D_NO_MPI $(DEFINES)
CXXFLAGS_MPI     = $(CXX_BASE) $(DEFINES)
CXXFLAGS_MPI    += -I$(METIS_INC_DIR) -I$(MPI_INC_DIR) 
CXXFLAGS_MPI    += -L$(METIS_LIB_DIR)

####### Output directory - these do nothing currently

OBJECTS_DIR   = ./obj
DESTDIR       = ./bin

####### Files

SOURCES = 	src/global.cpp \
		src/matrix.cpp \
		src/input.cpp \
		src/ele.cpp \
		src/polynomials.cpp \
		src/operators.cpp \
		src/geo.cpp \
		src/output.cpp \
		src/face.cpp \
		src/intFace.cpp \
		src/boundFace.cpp \
		src/mpiFace.cpp \
		src/flux.cpp \
		src/flurry.cpp \
		src/solver.cpp

OBJECTS = 	obj/global.o \
		obj/matrix.o \
		obj/input.o \
		obj/ele.o \
		obj/polynomials.o \
		obj/operators.o \
		obj/geo.o \
		obj/output.o \
		obj/face.o \
		obj/intFace.o \
		obj/boundFace.o \
		obj/mpiFace.o \
		obj/flux.o \
		obj/flurry.o \
		obj/solver.o

TARGET        = Flurry

####### Implicit rules

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(DESTDIR)/$(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS) $(DBG)

clean:
	cd obj && rm *.o && cd .. && rm bin/Flurry

.PHONY: debug
debug: CXXFLAGS=$(CXXFLAGS_DEBUG)
debug: LFLAGS+=-pg
debug: $(TARGET)

.PHONY: release
release: CXXFLAGS=$(CXXFLAGS_RELEASE)
release: $(TARGET)

.PHONY: openmp
openmp: CXXFLAGS=$(CXXFLAGS_OPENMP)
openmp: LIBS+= -fopenmp -lgomp
openmp: $(TARGET)

.PHONY: mpi
mpi: CXX=$(MPICXX)
mpi: LINK=$(MPILD)
mpi: CXXFLAGS=$(CXXFLAGS_MPI) $(CXX_RELEASE)
mpi: LIBS+= -lmetis
mpi: $(TARGET)

.PHONY: mpidebug
mpidebug: CXX=$(MPICXX)
mpidebug: LINK=$(MPILD)
mpidebug: CXXFLAGS=$(CXXFLAGS_MPI) $(CXX_DEBUG)
mpidebug: LIBS+= -lmetis
mpidebug: $(TARGET)

####### Compile

obj/global.o: src/global.cpp include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/global.o src/global.cpp

obj/matrix.o: src/matrix.cpp include/matrix.hpp \
		include/error.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/matrix.o src/matrix.cpp

obj/input.o: src/input.cpp include/input.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/input.o src/input.cpp

obj/ele.o: src/ele.cpp include/ele.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/flux.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ele.o src/ele.cpp

obj/polynomials.o: src/polynomials.cpp include/polynomials.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/polynomials.o src/polynomials.cpp

obj/operators.o: src/operators.cpp include/operators.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/face.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/operators.o src/operators.cpp

obj/geo.o: src/geo.cpp include/geo.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/geo.inl
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/geo.o src/geo.cpp

obj/output.o: src/output.cpp include/output.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/output.o src/output.cpp

obj/face.o: src/face.cpp include/face.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/ele.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/flux.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/face.o src/face.cpp

obj/intFace.o: src/intFace.cpp include/intFace.hpp  include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/intFace.o src/intFace.cpp

obj/boundFace.o: src/boundFace.cpp include/boundFace.hpp include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/boundFace.o src/boundFace.cpp

obj/mpiFace.o: src/mpiFace.cpp include/mpiFace.hpp include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/mpiFace.o src/mpiFace.cpp

obj/flux.o: src/flux.cpp include/flux.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/input.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/flux.o src/flux.cpp

obj/flurry.o: src/flurry.cpp include/flurry.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/input.hpp \
		include/geo.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/output.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/flurry.o src/flurry.cpp

obj/solver.o: src/solver.cpp include/solver.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/ele.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/solver.o src/solver.cpp
