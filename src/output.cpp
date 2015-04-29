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

#include <iomanip>

void writeData(solver *Solver, input *params, int iter)
{
  ofstream dataFile;

  char fileNameC[50];
  string fileName = params->dataFileName;
  //sprintf(fileNameC,"%s_%.09d.csv",&fileName[0],iter);
  sprintf(fileNameC,"%s.csv.%.09d",&fileName[0],iter);

  dataFile.open(fileNameC);

  // Vector of primitive variables
  vector<double> V;
  // Location of solution point
  point pt;

  // Eventually I'll get around to a real .vtk output (or Tecplot output...)
  // For now, let's output in a simple, Matlab-readable format: For each solution point:

  // Write header:
  // x  y  z(=0)  rho  [u  v  p]
  dataFile << "x,y,z,";
  if (params->equation == ADVECTION_DIFFUSION) {
    dataFile << "rho" << endl;
  }
  else if (params->equation == NAVIER_STOKES) {
    dataFile << "rho,u,v,p" << endl;
  }

  // Solution data
  for (auto& e:Solver->eles) {
    for (uint spt=0; spt<e.getNSpts(); spt++) {
      V = e.getPrimitives(spt);
      pt = e.getPosSpt(spt);

      for (uint dim=0; dim<e.getNDims(); dim++) {
        dataFile << pt[dim] << ",";
      }
      if (e.getNDims() == 2) dataFile << "0.0,"; // output a 0 for z [2D]

      for (uint i=0; i<e.getNFields()-1; i++) {
        dataFile << V[i] << ",";
      }
      dataFile << V[e.getNFields()-1] << endl;
    }
  }

  dataFile.close();
}


void writeResidual(solver *Solver, input *params, int iter)
{
  vector<double> res(params->nFields), resTmp(params->nFields);

  for (auto& e:Solver->eles) {
    resTmp = e.getResidual(params->resType);
    if(checkNaN(resTmp)) FatalError("NaN Encountered in Solution Residual!");

    for (int i=0; i<params->nFields; i++) {
      if (params->resType == 3) {
        // Infinity norm [max residual over all spts]
        res[i] = max(res[i],resTmp[i]);
      }else{
        res[i] += resTmp[i];
      }
    }
  }

  // If taking 2-norm, res is sum squared; take sqrt to complete
  if (params->resType == 2) {
    for (auto& i:res) i = sqrt(i);
  }

  int colW = 16;
  cout.precision(10);
  cout.setf(ios::fixed, ios::floatfield);
  if (iter==1 || (iter/params->monitor_res_freq)%40==0) {
    cout << endl;
    cout << setw(8) << left << "Iter";
    if (params->equation == ADVECTION_DIFFUSION) {
      cout << " Residual " << endl;
    }else if (params->equation == NAVIER_STOKES) {
      cout << setw(colW) << left << "rho";
      cout << setw(colW) << left << "rhoU";
      cout << setw(colW) << left << "rhoV";
      cout << setw(colW) << left << "rhoE";
    }
    cout << endl;
  }

  cout << setw(8) << left << iter;
  for (int i=0; i<params->nFields; i++) {
    cout << setw(colW) << left << res[i];
  }
  cout << endl;
}
