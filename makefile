#############################################################################
# Makefile for building Flurry
# Command: make debug
#          make release
#          make openmp
#############################################################################

#include configfiles/default.config
include configfiles/centos_intel.config

####### Compiler, tools and options

LIBS          = $(SUBLIBS)

TIOGA_INC   = ./lib/tioga/src

CXX_BASE    = -pipe -Wunused-parameter -Wuninitialized -std=c++11 -I./include -I$(TIOGA_INC) $(DEFINES)

CXX_BLAS = -I$(BLAS_INC_DIR) -L$(BLAS_LIB_DIR)
ifeq ($(strip $(BLAS_TYPE)),ATLAS)
	# --- ATLAS BLAS ---
  LD_BLAS = -L$(BLAS_LIB_DIR) -latlas -lcblas
else
ifeq ($(strip $(BLAS_TYPE)),MKL)
	# --- MKL BLAS ---
  ifeq ($(strip $(MPI)),YES)
	  # Serial BLAS [For use with MPI]
	  LD_BLAS = -I$(BLAS_INC_DIR) -L$(BLAS_LIB_DIR) $(BLAS_LIB_DIR)/libmkl_intel_lp64.a -Wl,--start-group $(BLAS_LIB_DIR)/libmkl_sequential.a $(BLAS_LIB_DIR)/libmkl_core.a -Wl,--end-group -L$(BLAS_COMPLIB_DIR) -lpthread -ldl -lm
  else
	  # Multithreaded BLAS [DO NOT USE WITH MPI]
	  LD_BLAS = -I$(BLAS_INC_DIR) -L$(BLAS_LIB_DIR) $(BLAS_LIB_DIR)/libmkl_intel_lp64.a -Wl,--start-group $(BLAS_LIB_DIR)/libmkl_intel_thread.a $(BLAS_LIB_DIR)/libmkl_core.a -Wl,--end-group -L$(BLAS_COMPLIB_DIR) -liomp5 -lpthread -ldl -lm
  endif
	CXX_BLAS += -D_MKL_BLAS
else
	# --- Other BLAS ---
  LD_BLAS = -L$(BLAS_LIB_DIR) -lcblas
endif
endif

CXX_BASE += $(CXX_BLAS)
LIBS += $(LD_BLAS)

CXXFLAGS = $(CXX_BASE)
ifeq ($(strip $(DEBUG_LEVEL)),0)
  CXXFLAGS += -Ofast -fno-finite-math-only
  FFLAGS = -Ofast
else ifeq ($(strip $(DEBUG_LEVEL)),1)
  CXXFLAGS += -g -O2
  FFLAGS = -Ofast
else ifeq ($(strip $(DEBUG_LEVEL)),2)
  CXXFLAGS += -g -pg -O0 -rdynamic -fno-omit-frame-pointer #-fsanitize=address 
  FFLAGS = -g -O0 -rdynamic -fno-omit-frame-pointer #-fsanitize=address
endif

ifeq ($(strip $(ENABLE_DEBUG)),1)
  CXXFLAGS += -D_DEBUG
endif

ifeq ($(strip $(OPENMP)),YES)
  CXXFLAGS += -fopenmp -D_OMP
  LFLAGS += -fopenmp
else
  CXXFLAGS += -Wno-unknown-pragmas
endif

ifeq ($(strip $(MPI)),YES)
  CXXFLAGS += -Wno-literal-suffix -I$(METIS_INC_DIR) -I$(MPI_INC_DIR) -L$(METIS_LIB_DIR)
  ifeq ($(MPI_DEBUG),YES)
    CXXFLAGS += -D_MPI_DEBUG
  endif
  LIBS += -L$(METIS_LIB_DIR) -lmetis
  CXX = $(MPICXX)
  LINK = $(MPILD)
else
  CXXFLAGS += -D_NO_MPI
endif

####### Output directory - these do nothing currently

OBJECTS_DIR   = ./obj
DESTDIR       = ./bin

####### Files

OBJECTS = 	obj/global.o \
		obj/funcs.o \
		obj/points.o \
		obj/matrix.o \
		obj/input.o \
		obj/ele.o \
		obj/polynomials.o \
		obj/operators.o \
		obj/geo.o \
		obj/geo_overset.o \
		obj/output.o \
		obj/face.o \
		obj/intFace.o \
		obj/boundFace.o \
		obj/mpiFace.o \
		obj/overFace.o \
		obj/flux.o \
		obj/flurry.o \
		obj/solver.o \
		obj/solver_overset.o \
		obj/multigrid.o \
		obj/superMesh.o \
		obj/overComm.o

ifeq ($(strip $(MPI)),YES)
OBJECTS+= obj/ADT.o \
		obj/MeshBlock.o \
		obj/parallelComm.o \
		obj/utils.o \
		obj/tioga.o \
		obj/kaiser.o \
		obj/median.o \
		obj/cellVolume.o
endif

TARGET        = Flurry

####### Implicit rules

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.f90.o:
	$(F90) -c $(FFLAGS) $(INCPATH) -o "$@" "$<"
	
.F90.o:
	$(F90) -c $(FFLAGS) $(INCPATH) -o "$@" "$<"

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(DESTDIR)/$(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS) $(DBG)

####### Build rules

clean:
	cd obj && rm -f *.o && cd .. && rm -f bin/Flurry

####### Compile

obj/global.o: src/global.cpp include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/global.o src/global.cpp

obj/funcs.o: src/funcs.cpp include/funcs.hpp include/global.hpp include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/funcs.o src/funcs.cpp

obj/points.o: src/points.cpp include/points.hpp include/funcs.hpp include/global.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/points.o src/points.cpp

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
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/geo.o src/geo.cpp

obj/geo_overset.o: src/geo_overset.cpp include/geo.hpp \
	include/ele.hpp \
	include/face.hpp \
	include/overFace.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/geo_overset.o src/geo_overset.cpp

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

obj/overFace.o: src/overFace.cpp include/overFace.hpp include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/overFace.o src/overFace.cpp

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
		include/overComm.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/solver.o src/solver.cpp

obj/solver_overset.o: src/solver_overset.cpp include/solver.hpp \
  		include/global.hpp \
		include/input.hpp \
		include/geo.hpp \
		include/matrix.hpp \
		include/overComm.hpp \
		include/overFace.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/solver_overset.o src/solver_overset.cpp

obj/multigrid.o: src/multigrid.cpp include/multigrid.hpp \
	include/input.hpp \
	include/solver.hpp \
	include/ele.hpp \
	include/geo.hpp \
	include/operators.hpp \
	include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/multigrid.o src/multigrid.cpp

obj/superMesh.o: src/superMesh.cpp include/superMesh.hpp \
	include/global.hpp \
	include/matrix.hpp \
	include/geo.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/superMesh.o src/superMesh.cpp
	
obj/overComm.o: src/overComm.cpp include/overComm.hpp \
		include/global.hpp \
		include/funcs.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/matrix.hpp \
		include/operators.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/overComm.o src/overComm.cpp


obj/ADT.o: lib/tioga/src/ADT.C lib/tioga/src/ADT.h \
  	lib/tioga/src/codetypes.h \
	lib/tioga/src/utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ADT.o lib/tioga/src/ADT.C

obj/utils.o: lib/tioga/src/utils.c lib/tioga/src/utils.h 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/utils.o lib/tioga/src/utils.c
	
obj/parallelComm.o: lib/tioga/src/parallelComm.C lib/tioga/src/parallelComm.h \
  	lib/tioga/src/codetypes.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/parallelComm.o lib/tioga/src/parallelComm.C
	
obj/MeshBlock.o: lib/tioga/src/MeshBlock.C lib/tioga/src/MeshBlock.h \
  	lib/tioga/src/codetypes.h \
	lib/tioga/src/utils.h \
	lib/tioga/src/ADT.h \
	lib/tioga/src/parallelComm.h \
	include/solver.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/MeshBlock.o lib/tioga/src/MeshBlock.C

obj/tioga.o: lib/tioga/src/tioga.C lib/tioga/src/tioga.h \
  	lib/tioga/src/codetypes.h \
	lib/tioga/src/utils.h \
	lib/tioga/src/parallelComm.h \
	lib/tioga/src/MeshBlock.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/tioga.o lib/tioga/src/tioga.C

obj/kaiser.o: lib/tioga/src/kaiser.f 
	$(F90) -c $(FFLAGS) $(INCPATH) -o obj/kaiser.o lib/tioga/src/kaiser.f
	
obj/cellVolume.o: lib/tioga/src/cellVolume.f90
	$(F90) -c $(FFLAGS) $(INCPATH) -o obj/cellVolume.o lib/tioga/src/cellVolume.f90

obj/median.o: lib/tioga/src/median.F90 
	$(F90) -c $(FFLAGS) $(INCPATH) -o obj/median.o lib/tioga/src/median.F90
