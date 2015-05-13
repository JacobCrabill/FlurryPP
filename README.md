FlurryPP
========

A 2D Flux Reconstruction Code in C++

Written by Jacob Crabill
Aerospace Computing Lab, Stanford University

Current Capabilities
====================

The code is currently capable of running scalar advection or Euler (inviscid Navier-Stokes) cases on unstructured mixed grids of quadrilaterals and triangles in the Gmsh format.
Moving grids are supported by the solver, but there are not yet any grid-motion functions implemented beyond a standard test case.

Background / Goals of the Project
=================================

This is not intended to be a high-performance CFD tool, but rather a learning experience and testbed for new ideas and methods.  The idea is to have a simple, easy to understand codebase which has been created from the ground up to be easy to modify in the future.  As a developer of HiFiLES, I have great respect for the high-performance aspect of the code, and the thought that went into optimizing its performance for a specific application.  However, I've also learned the hard way that it's rather difficult to make an established code do things that it wasn't designed to do.

Hence, Flurry was born.  To start with, my goal is to create a plain 2D Euler and Navier-Stokes solver on quadrilateral and/or triangular elements, mostly as a learning experience.  As more long-term goals, however, I hope to implement various mesh adaptation and deformation methods, p-adaptation, h- and p-multigrid, etc., and as such the code will be structured in such a way to make that possible some day.

In the meantime, if you actually took the time to read this, then you should really check out the high-performance GPU-capable Flux Reconstruction code HiFiLES: http://github.com/HiFiLES/HiFiLES-solver
This open-source code is under development by the Aerospace Computing Lab at Stanford, with new features in the works all the time (it is a research code, after all).


Quick Start
===========

Compilation Instructions
-------------------------

To compile Flurry, you can either use QT Creator (https://www.qt.io/download-open-source/), which is an excellent C++ IDE that I use for development, or you can use the provided Makefile.flurry to compile using GNU make.  For the make option, just open a terminal and run the following:

`make`

Optionally, you can specify the type of build as either *debug* or *release*:

`make release`

where *release* turns on full optimization, and *debug* removes all optimization and adds flags for both debugging and profiling.  
The code also (optionally) utilizes OpenMP to take advantage of easy parallelization on desktop computers; speedup of ~3x has been seen on 4-core processors.  To enable OpenMP when compiling, just do:

`make openmp`


Test Cases
-------------------------

Several basic test cases have been created to verify the functionality of the code.  The first test case is for the advection equation (in 'tests/advection'), and is simply the advection of a Gaussian bump in a periodic domain.
The other two test cases are for the inviscid Navier-Stokes (Euler) equations.  
One is for supersonic flow over a wedge, and the other is for subsonic flow over a circular cylinder.  
Note that Flurry does not currently have any shock-capturing methods implemented, so although this test case works, general transonic and supersonic cases should be approached with caution.
The cylinder test case uses a very coarse mesh, and is intended purely for the purpose of testing the functionality of the code on arbitrary unstructured quad meshes from Gmsh, and demonstrating the method for applying boundary conditions to Gmsh meshes.


Post-Processing
-------------------------

Flurry currently has two options for viewing simulation data: ParaView .vtu files, and a super-simple .csv file.  In both cases, the output values are the primitive variables, not the conservative variables.  

ParaView is a free, cross-platform visualization tool that works quite well for visualizing CFD data; you can get it from http://www.paraview.org/.

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
- Solver
  + Applies the various FR operations to a solution (set of eles, faces, operators, and geometry)


Potential Classes (Food for thought)
------------------------------------
- Solution
  + Wrapper for all data needed for a solution; for h-multigrid or overset, can have several, with methods to transfer data between the two?
  + Alternatively - have a "solution" class/struct which just holds the global solution, flux, common flux vectors, and then have the eles (and inters?) classes link to this global array somehow (hopefully, something a little more clean/elegant that HiFiLES's "inters" class; sure, pointers work, but having to dereference a pointer on every single variable gets annoying)


Future Classes (?)
------------------
- Overset
  + given 2 grids, find the overlapping region and handle data transfer between the two solutions
