/*!
 * \file output.cpp
 * \brief Functions for solution restart & visualization data output
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
#include "../include/output.hpp"

void writeData(solver *Solver, input *params, int iter)
{
  ofstream dataFile;

  char fileNameC[50];
  string fileName = params->dataFileName;
  sprintf(fileNameC,"%s_%.09d.vtk",&fileName[0],iter);

  dataFile.open(fileNameC);

  // Vector of primitive variables
  vector<double> V;
  // Location of solution point
  point pt;

  // Eventually I'll get around to a real .vtk output (or Tecplot output...)
  // For now, let's output in a simple, Matlab-readable format: For each solution point:
  // x  y  rho  u  v  p
  for (auto& e:Solver->eles) {
    for (uint spt=0; spt<e.getNSpts(); spt++) {
      V = e.getPrimitives(spt);
      pt = e.getPosSpt(spt);
      for (int dim=0; dim<e.getNDims(); dim++) {
        dataFile << pt[dim] << ",";
      }
      for (int i=0; i<e.getNFields()-1; i++) {
        dataFile << V[i] << ",";
      }
      dataFile << V[e.getNFields()-1] << endl;
    }
  }

  dataFile.close();
}
