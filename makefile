#############################################################################
# Makefile for building Flurry
#############################################################################

include configfiles/default.config
#include configfiles/centos_intel.config

####### Compiler, tools and options

SWIG = swig -c++ -python
LIBS = $(SUBLIBS)
TIOGA_INC = ./lib/tioga/src

CXX_BASE = -pipe -Wunused-parameter -Wuninitialized -std=c++11 
INCS = -I./include -I$(TIOGA_INC) $(DEFINES)

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
ifeq ($(strip $(BLAS_TYPE)),OPEN)
  LD_BLAS = -L$(BLAS_LIB_DIR) -lopenblas
else
	# --- Other BLAS ---
  LD_BLAS = -L$(BLAS_LIB_DIR) -lcblas
endif
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
  FLAGS += -D_DEBUG
endif

ifeq ($(strip $(OPENMP)),YES)
  CXXFLAGS += -fopenmp
  FLAGS += -D_OMP
  LFLAGS += -fopenmp
else
  CXXFLAGS += -Wno-unknown-pragmas
endif

ifeq ($(strip $(MPI)),YES)
  CXXFLAGS += -Wno-literal-suffix -L$(METIS_LIB_DIR)
	INCS += -I$(METIS_INC_DIR) -I$(MPI_INC_DIR) 
  ifeq ($(MPI_DEBUG),YES)
    FLAGS += -D_MPI_DEBUG
  endif
  LIBS += -L$(METIS_LIB_DIR) -lmetis
  CXX = $(MPICXX)
  LINK = $(MPILD)
else
  FLAGS += -D_NO_MPI
endif

####### Output directory - these do nothing currently

OBJDIR        = $(CURDIR)/obj
BINDIR       = $(CURDIR)/bin
SWIGDIR       = $(CURDIR)/swig

####### Files

OBJECTS = 	$(OBJDIR)/global.o \
		$(OBJDIR)/funcs.o \
		$(OBJDIR)/points.o \
		$(OBJDIR)/matrix.o \
		$(OBJDIR)/input.o \
		$(OBJDIR)/ele.o \
		$(OBJDIR)/polynomials.o \
		$(OBJDIR)/operators.o \
		$(OBJDIR)/geo.o \
		$(OBJDIR)/geo_overset.o \
		$(OBJDIR)/output.o \
		$(OBJDIR)/face.o \
		$(OBJDIR)/intFace.o \
		$(OBJDIR)/boundFace.o \
		$(OBJDIR)/mpiFace.o \
		$(OBJDIR)/overFace.o \
		$(OBJDIR)/flux.o \
		$(OBJDIR)/flurry.o \
		$(OBJDIR)/solver.o \
		$(OBJDIR)/solver_overset.o \
		$(OBJDIR)/multigrid.o \
		$(OBJDIR)/superMesh.o \
		$(OBJDIR)/overComm.o

ifeq ($(strip $(MPI)),YES)
OBJECTS+= $(OBJDIR)/ADT.o
#		$(OBJDIR)/tioga.o 
#		$(OBJDIR)/MeshBlock.o \
#		$(OBJDIR)/parallelComm.o \
#		$(OBJDIR)/utils.o \
#		$(OBJDIR)/kaiser.o \
#		$(OBJDIR)/median.o \
#		$(OBJDIR)/cellVolume.o
endif

WRAP_OBJS = $(OBJECTS) $(OBJDIR)/flurry_interface.o
SWIG_OBJS = $(OBJECTS) $(OBJDIR)/flurry_interface.o $(SWIGDIR)/flurry_wrap.o
SWIG_INCS = -I/usr/include/python2.7 -I/usr/local/lib/python2.7/dist-packages/mpi4py/include/
SWIG_LIBS =
INCS += $(SWIG_INCS)

TARGET        = Flurry

####### Implicit rules
$(TARGET): $(OBJECTS)

	$(LINK) $(INCS) $(LFLAGS) -o $(BINDIR)/$(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS) $(DBG)

.PHONY: swig
swig: CXXFLAGS += -fPIC -D_BUILD_LIB $(INCS)
swig: FFLAGS += -fPIC
swig: $(SWIG_OBJS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INCS) $(SWIG_INCS) -shared -o $(SWIGDIR)/_flurry.so $(SWIG_OBJS) $(LIBS) $(SWIG_LIBS)

.PHONY: lib
lib: CXXFLAGS += -fPIC -D_BUILD_LIB
lib: FFLAGS += -fPIC
lib: $(WRAP_OBJS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INCS) $(SWIG_INCS) -shared -o $(BINDIR)/libflurry.so $(WRAP_OBJS) $(LIBS) $(SWIG_LIBS)

.PHONY: test
test: INCS += -I~/tioga/src/
test: lib
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCS) $(SWIGDIR)/testFlurry.cpp $(BINDIR)/libflurry.so -L$(SWIGDIR)/lib -ltioga -Wl,-rpath=$(SWIGDIR)/lib/ -o $(SWIGDIR)/testFlurry
	cp $(BINDIR)/libflurry.so $(SWIGDIR)/lib/

####### Implicit Rules

$(OBJDIR)/%.o: src/%.cpp 
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INCS) $(SWIG_INCS) -c -o $@ $<

$(OBJDIR)/%.o: $(SWIGDIR)/%.cpp 
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INCS) $(SWIG_INCS) -c -o $@ $<

$(OBJDIR)/%.o: lib/tioga/src/%.C 
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INCS) $(SWIG_INCS) -c -o $@ $<

$(OBJDIR)/%.o: lib/tioga/src/%.c 
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INCS) $(SWIG_INCS) -c -o $@ $<

$(OBJDIR)/%.o: lib/tioga/src/%.f	
	$(F90) -c $(INCS) $(FFLAGS) -o $@ $< 

$(OBJDIR)/%.o: lib/tioga/src/%.f90
	$(F90) -c $(INCS) $(FFLAGS) -o $@ $< 

$(OBJDIR)/%.o: lib/tioga/src/%.F90
	$(F90) -c $(INCS) $(FFLAGS) -o $@ $< 
	
$(SWIGDIR)/%_wrap.cpp: $(SWIGDIR)/%.i
	$(SWIG) $(INCS) $(FLAGS) $(SWIG_INCS) -o $@ $<

####### Build rules

clean:
	cd $(OBJDIR) && rm -f *.o && cd .. && rm -f bin/Flurry
	rm -f lib/tioga/src/*.o


