%module flurry

// -----------------------------------------------------------------------------
// Header files required by any of the following C++ code
// -----------------------------------------------------------------------------
%header
%{
#include "flurry_interface.hpp"
%}

// -----------------------------------------------------------------------------
// Header files and other declarations to be parsed as SWIG input
// -----------------------------------------------------------------------------

#ifndef _NO_MPI
// MPI SWIG interface file & MPI_Comm to Python Comm typemap
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm,MPI_Comm);
#endif

%include "input.hpp"
%include "flurry.hpp"
%include "flurry_interface.hpp"

// <-- Additional C++ declations [anything that would normally go in a header]

// -----------------------------------------------------------------------------
// Additional functions which have been declared, but not defined (including
// definition in other source files which will be linked in later)
// -----------------------------------------------------------------------------

%inline
%{
// <-- Additional C++ definitions [anything that would normally go in a .cpp]
%}

// -----------------------------------------------------------------------------
// Additional Python functions to add to module
// [can use any functions/variables declared above]
// -----------------------------------------------------------------------------

%pythoncode
%{
# Python functions here
%}


