/*!
 * \file output.hpp
 * \brief Header file for restart & visualization data output functions
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include "global.hpp"
#include "solver.hpp"
#include "ele.hpp"
#include "geo.hpp"

/*! Write solution to file (of type params->plotType) */
void writeData(solver *Solver, input *params);

/*! Write solution data to a CSV file. */
void writeCSV(solver *Solver, input *params);

/*! Write solution data to a Paraview .vtu file. */
void writeParaview(solver *Solver, input *params);

/*! Compute the residual and print to both the terminal and history file. */
void writeResidual(solver *Solver, input *params);

/*! Compute and display all error norms */
void writeAllError(solver *Solver, input *params);

/*! Compute the conservation, L1, or L2 solution error (for certain test cases) and print to screen. */
void writeError(solver *Solver, input *params);

/*! Write a Tecplot mesh file compatible with TIOGA's testTioga FORTRAN interface */
void writeMeshTecplot(solver* Solver, input* params);
