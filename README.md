FlurryPP
========

A 2D & 3D Flux Reconstruction Code in C++

Written by Jacob Crabill

Aerospace Computing Lab, Stanford University

Current Capabilities
====================

The code is currently capable of running scalar advection/diffusion or Euler/Navier-Stokes cases on unstructured mixed grids of quadrilaterals and triangles (2D) or hexahedrons (3D) in the Gmsh format.  Both linear and quadratic elements are supported in 3D, and up to 10th-order quadrilaterals in 2D, allowing for accurate representations of curved boundaries.

Local/Global CFL-based time stepping is available for ease (and safety) of use, along with both Forward Euler and RK45 time-stepping.

Convergence accleration is available with the use of P- or HP-multigrid.  The use of H multigrid is only in addition to P multigrid down to order 0, and uses a somewhat novel grid-generation technique in which the solver accepts the coarsest grid level to run on, and recursively refines the mesh to generate the intermediate levels up to a final refined grid.  See the tests/euler/channel test case for examples of using both MG methods.

Shock capturing has been implemented, but is still under development and is not fully tested yet. Additionally, a highly robust stabilization procedure invented by Chi-Wang Shu and further developed by Yu Lv is available; however, its usage tends to disrupt convergence of steady-state problems.

Moving grids are supported by the solver, but there are not yet any grid-motion functions implemented beyond several simple test-case functions.

Lastly, overset grids in 2D and 3D are supported by using the "artificial boundary" method, with Jay Sitaraman's TIOGA library being used for hole-blanking whenever solid bodies are embedded inside a mesh.  Moving overset grids are also supported now for 2D and 3D.

Background / Goals of the Project
=================================

This is not intended to be a high-performance CFD tool, but rather a learning experience and testbed for new ideas and methods.  The idea is to have a simple, easy to understand codebase which has been created from the ground up to be easy to modify in the future.  As a developer of HiFiLES, I have great respect for the high-performance aspect of the code, and the thought that went into optimizing its performance for a specific application.  However, I've also learned the hard way that it's rather difficult to make an established code do things that it wasn't designed to do.

Hence, Flurry was born.  To start with, my goal is to create a plain 2D/3D Euler and Navier-Stokes solver on quadrilateral and/or triangular elements, mostly as a learning experience.  As more long-term goals, however, I hope to implement various mesh adaptation and deformation methods, p-adaptation, h- and p-multigrid, etc., and as such the code will be structured in such a way to make that possible some day.

In the meantime, if you actually took the time to read this, then you should really check out the high-performance GPU-capable Flux Reconstruction code HiFiLES: http://github.com/HiFiLES/HiFiLES-solver
This open-source code is under development by the Aerospace Computing Lab at Stanford, with new features in the works all the time (it is a research code, after all).


Quick Start
===========

Compilation Instructions
-------------------------

To compile Flurry, you can use the provided makefile to compile using GNU make.  First modify 
the config file (configfiles/default.config) or create your own with your desired compilation 
options, then from FlurryPP directory, type `make`.</p>

See the default config file for details on all compilation options. Note that compiling with MPI 
requires several external libraries and header files (metis.h/metis.a, mpi.h), the location of 
which must be specified in the makefile.
The code is known to work with OpenMPI >= 1.6.5.
Additionally, you should have some form of BLAS installed, and ensure that the location of cblas.h 
is properly specified in the config file.


Test Cases
-------------------------

Several basic test cases have been created to verify the functionality of the code.  The first test case is for the advection equation (in 'tests/advection'), and is simply the advection of a Gaussian bump in a periodic domain.

The other two test cases are for the inviscid Navier-Stokes (Euler) equations. One is for supersonic flow over a wedge, another is for subsonic flow over a circular cylinder, and another is for inviscid channel flow over a Gaussian bump..  
Note that although Flurry does have a shock-capturing method implemented, it is still in the developmental phase, so general transonic and supersonic cases should be approached with caution.
The cylinder test case uses a very coarse mesh, and is intended purely for the purpose of testing the functionality of the code on arbitrary unstructured quad meshes from Gmsh, and demonstrating the method for applying boundary conditions to Gmsh meshes.

Lastly, there are several test cases available for overset grids in both 2D and 3D to see the available functionality.


Post-Processing
-------------------------

Flurry currently has two options for viewing simulation data: ParaView .vtu files, and a super-simple .csv file.  In both cases, the output values are the primitive variables, not the conservative variables.

ParaView is a free, cross-platform visualization tool that works quite well for visualizing CFD data, both 2D and 3D; you can get the latest version from http://www.paraview.org/.

The CSV output method, on the other hand, simply outputs the x,y,z coordinates of the solution points in each element, along with the solution vector at each point. This can either be plotted using the provided Matlab script.
If you write any similar scripts for plotting in other languages (e.g. Python or Julia), please let me know so I can add them here!


Code Structure
==============

Files
-----
- Flurry
  + Driver
- Input
  + Read input file, set run parameters
- Output
  + Write data to file for restarting, visualization, etc.
- Geometry
  + Read mesh file / generate mesh; setup eles & faces
- Flux
  + Routines to calculate the inviscid & viscous fluxes
- Polynomials
  + Lagrange, Legendre, and other useful polynomials and functions


Basic Classes
--------------
- Operators
  + Pre-computes and stores matrices for interpolation, extrapolation, etc. of polynomial bases for all elements in current solution
- Ele
  + Each 'ele' is a single element in the mesh, and stores its solution, Jacobian, shape nodes, and global face IDs of its faces
- Face
  + Each face stores the ID of the L/R cells & local face ID of each
  + Calculates interface flux or boundary flux
  + Passes the interface flux difference (common minus discontinuous) back to the L/R elements
- Bound
  + Each 'bound' is just like a face, except that the 'right' element is replaced by a 'ghost' right state from boundary conditions
  + The boundary conditions are applied to the right state in such a way that a central flux between the left and right states produces the necessary common flux
- MPI Face
  + Nearly identical to an internal face, except that communication of right-state data is handled through non-blocak MPI sends/recieves.
  + Note that, as opposed to internal faces, MPI faces are duplicated across processor boundaries
- Overset Face
  + Similar in concept to the MPI face, all overset faces grab their right state from the solver to which they belong, after the solver has performed the interpolation from the other grid(s).
- OverComm (Overset Communicator)
  + Handles all data communication related to overset grids, and implements some 2D hole-cutting algorithms which do not exist in Tioga.
- SuperMesh
  + Creates a local supermesh of a target element from one or more donor elements for use with Galerkin projection (see Farrell and Maddison, 2010).
- Solver
  + Applies the various FR operations to a solution (set of eles, faces, operators, and geometry)
- MultiGrid
  + Sets up the additional solvers needed to run P- or HP-multigrid, and performs a V-cycle multigrid update.


Flurry++ is distributed under the [GNU General Public License Version 3](http://www.gnu.org/licenses/gpl-3.0.txt).

![GPL3](http://www.gnu.org/graphics/gplv3-127x51.png)
