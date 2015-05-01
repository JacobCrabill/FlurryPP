/*!
 * \file output.hpp
 * \brief Header file for restart & visualization data output functions
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */
#pragma once

#include "global.hpp"
#include "solver.hpp"
#include "ele.hpp"
#include "geo.hpp"

/*! Write solution to file (of type params->plot_type) */
void writeData(solver *Solver, input *params);

/*! Write solution data to a CSV file. */
void writeCSV(solver *Solver, input *params);

/*! Write solution data to a Paraview .vtu file. */
void writeParaview(solver *Solver, input *params);
void writeParaviewBinary(solver *Solver, input *params);

/*! Compute the residual and print to the screen. */
void writeResidual(solver *Solver, input *params);
