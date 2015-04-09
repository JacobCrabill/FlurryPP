FlurryPP
========

A 2D Flux Reconstruction Code in C++

Written by Jacob Crabill
Aerospace Computing Lab, Stanford University


Background / Goals of the Project
=================================

This is not intended to be a high-performance CFD tool, but rather a learning experience and testbed for new ideas and methods.  The idea is to have a simple, easy to understand codebase which has been created from the ground up to be easy to modify in the future.  As a developer of HiFiLES, I have great respect for the high-performance aspect of the code, and the thought that went into optimizing its performance for a specific application.  However, I've also learned the hard way that it's rather difficult to make an established code do things that it wasn't designed to do.

Hence, Flurry was born.  To start with, my goal is to create a plain 2D Euler and Navier-Stokes solver on quadrilateral and/or triangular elements, mostly as a learning experience.  As more long-term goals, however, I hope to implement various mesh adaptation and deformation methods, p-adaptation, h- and p-multigrid, etc., and as such the code will be structured in such a way to make that possible some day.

In the meantime, if you actually took the time to read this, then you should really check out the high-performance GPU-capable Flux Reconstruction code HiFiLES: http://github.com/HiFiLES/HiFiLES-solver
This open-source code is under development by the Aerospace Computing Lab at Stanford, with new features in the works all the time (it is a research code, after all).


Quick Start - Compilation Instructions
=======================================

To compile Flurry, you can either use QT Creator (https://www.qt.io/download-open-source/), which is an excellent C++ IDE that I use for development, or you can use the provided Makefile.flurry to compile using GNU make.  For the make option, just open a terminal and run the following:
`make -f Makefile.flurry`
Optionally, you can specify the type of build as either *debug* or *release*:
`make -f Makefile.flurry CODE='release'`
where *release* turns on full optimization, and *debug* removes all optimization and adds flags for both debugging and profiling.


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


Basic Classes
--------------
- Operators
  + Pre-computes and stores matrices for interpolation, extrapolation, etc. of polynomial bases for all elements in current solution
- Ele
  + Each 'ele' is a single element in the mesh, and stores its solution, Jacobian, shape nodes (pointer to global vertices?), and pointer(?) to global face IDs
- Face
  + Each face stores the ID of the L/R cells & local face ID of each
  + Calculates interface flux or boundary flux
- Solver (grid? mesh? solution?)
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
