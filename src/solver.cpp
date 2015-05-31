/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
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

#include "../include/solver.hpp"

#include <sstream>
#include <omp.h>
#include <stdexcept>

class intFace;
class boundFace;

#include "../include/input.hpp"
#include "../include/geo.hpp"
#include "../include/intFace.hpp"
#include "../include/boundFace.hpp"
#include "../include/output.hpp"

solver::solver()
{
}

void solver::setup(input *params, geo *Geo)
{
  this->params = params;
  this->Geo = Geo;

  params->time = 0.;

  /* Setup the FR elements & faces which will be computed on */
  Geo->setupElesFaces(eles,faces,mpiFaces);

  if (params->restart) {
    readRestartFile();
  }
  else {
    setupElesFaces();
  }

  /* Setup the FR operators for computation */
  setupOperators();

  /* Additional Setup */

  // Time advancement setup
  switch (params->timeType) {
    case 0:
      nRKSteps = 1;
      RKb = {1};
      break;
    case 4:
      nRKSteps = 4;
      RKa = {.5, .5, 1.};
      RKb = {1./6., 1./3., 1./3., 1./6.};
      break;
    default:
      FatalError("Time-Stepping type not supported.");
  }

#ifndef _NO_MPI
  finishMpiSetup();
#endif
}


double solver::runSim(const vector<double> &X, int field)
{
  params->Kp = X[0];
  params->Kd = X[1];
  params->Ki = X[2];

  cout << endl;
  cout << setw(5) << left << "Kp: " << setw(8) << setprecision(4) << params->Kp << setw(2) << " ";
  cout << setw(5) << left << "Kd: " << setw(8) << setprecision(4) << params->Kd << setw(2) << " ";
  cout << setw(5) << left << "Ki: " << setw(8) << setprecision(4) << params->Ki << setw(2) << " ";
  cout << endl;

  params->hist << params->Kp << "," << params->Kd << "," << params->Ki << ",";

  simTimer runTime;
  runTime.startTimer();

  /* Apply the initial condition */
  resetBoundFaces();
  initializeSolution();

  params->dataFileName = params->outputPrefix + "_" + to_string(params->evals);
  double err;

  try {
    for (params->iter = params->initIter+1; params->iter < params->iterMax; params->iter++) {
      update();

      if ((params->iter)%params->monitor_res_freq == 0 || params->iter==1) writeResidual(this,params);
      if ((params->iter)%params->plot_freq == 0) writeData(this,params);

      if (params->iter%50==0) {
        double maxDU = 0;
#pragma omp parallel for reduction(max:maxDU)
        for (uint i=0; i<faces.size(); i++) {
          maxDU = std::max(maxDU,faces[i]->maxDU);
        }
        cout << "MaxDU = " << maxDU << endl;
      }
    }
    params->evals++;

    err = calcNormError(field);
    writeData(this,params);
  }
  catch (std::exception e) {
    cout << "NaN Encountered.  Setting function value to INFINITY." << endl;

    runTime.stopTimer();
    runTime.showTime();

    params->evals++;

    err = (double)INFINITY;
  }

  cout << endl;
  cout << "Final L2 Error: " << setprecision(8) << err << endl;
  cout << endl;

  params->hist << err << endl;

  return err;
}

void solver::resetBoundFaces(void)
{
  for (int i = Geo->nFaces; i<Geo->nFaces+Geo->nBndFaces; i++) {
    faces[i]->setupRightState();
  }
}

void solver::update(void)
{
  // For RK time-stepping, store the starting solution values
  if (nRKSteps>1)
    copyUspts_U0();

  /* Intermediate residuals for Runge-Kutta time integration */

  for (int step=0; step<nRKSteps-1; step++) {

    if (step == 1)
      params->rkTime = params->time;
    else
      params->rkTime = params->time + RKa[step-1]*params->dt;

    moveMesh(step);

    calcResidual(step);

    timeStepA(step);

  }

  /* Final Runge-Kutta time advancement step */

  if (nRKSteps == 1)
    params->rkTime = params->time;
  else
    params->rkTime = params->time + params->dt;

  moveMesh(nRKSteps-1);

  calcResidual(nRKSteps-1);

  // Reset solution to initial-stage values
  if (nRKSteps>1)
    copyU0_Uspts();

  for (int step=0; step<nRKSteps; step++) {
    timeStepB(step);
  }

  params->time += params->dt;
}


void solver::calcResidual(int step)
{
  extrapolateU();

#ifndef _NO_MPI
  doCommunication();
#endif

  if (params->viscous || params->motion) {

    calcGradU_spts();

  }

  calcInviscidFlux_spts();

  extrapolateNormalFlux();

  calcInviscidFlux_faces();

#ifndef _NO_MPI
  calcInviscidFlux_mpi();
#endif

  if (params->viscous) {

    correctU();

    extrapolateGradU();

    calcViscousFlux_faces();

#ifndef _NO_MPI
    calcViscousFlux_mpi();
#endif

  }

  calcFluxDivergence(step);

  correctDivFlux(step);
}


void solver::timeStepA(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].timeStepA(step,RKa[step]);
  }
}


void solver::timeStepB(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].timeStepB(step,RKb[step]);
  }
}

void solver::copyUspts_U0(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].copyUspts_U0();
  }
}

void solver::copyU0_Uspts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].copyU0_Uspts();
  }
}

void solver::extrapolateU(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applySptsFpts(eles[i].U_spts,eles[i].U_fpts);
  }
}

void solver::extrapolateUMpts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applySptsMpts(eles[i].U_spts,eles[i].U_mpts);
  }
}

void solver::extrapolateSMpts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applySptsMpts(eles[i].S_spts,eles[i].S_mpts);
  }
}

void solver::extrapolateSFpts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applySptsFpts(eles[i].S_spts,eles[i].S_fpts);
  }
}

void solver::calcInviscidFlux_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].calcInviscidFlux_spts();
  }
}

void solver::doCommunication()
{
#pragma omp parallel for
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->communicate();
  }
}

void solver::calcInviscidFlux_faces()
{
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->calcInviscidFlux();
  }
}

void solver::calcInviscidFlux_mpi()
{
#pragma omp parallel for
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->calcInviscidFlux();
  }
}

void solver::calcViscousFlux_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].calcInviscidFlux_spts();
  }
}

void solver::calcViscousFlux_faces()
{
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->calcViscousFlux();
  }
}


void solver::calcViscousFlux_mpi()
{
#pragma omp parallel for
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->calcViscousFlux();
  }
}


void solver::calcGradF_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyGradFSpts(eles[i].F_spts,eles[i].dF_spts);
  }
}

void solver::transformGradF_spts(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].transformGradF_spts(step);
    //opers[eles[i].eType][eles[i].order].applyTransformGradFSpts(eles[i].dF_spts,eles[i].JGinv_spts,eles[i].gridVel_spts);
  }
}

void solver::calcFluxDivergence(int step)
{
  if (params->motion) {

    /* Use non-conservation-form chain-rule transformation
     * (See AIAA paper 2013-0998 by Liang, Miyaji and Zhang) */

    calcGradF_spts();

    transformGradF_spts(step);

  }else{

    /* Standard conservative form */
    calcDivF_spts(step);

  }
}

void solver::calcDivF_spts(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyDivFSpts(eles[i].F_spts,eles[i].divF_spts[step]);
  }
}

void solver::extrapolateNormalFlux(void)
{
  if (params->motion) {
    /* Extrapolate physical normal flux */
#pragma omp parallel for
    for (uint i=0; i<eles.size(); i++) {
      opers[eles[i].eType][eles[i].order].applyExtrapolateFn(eles[i].F_spts,eles[i].norm_fpts,eles[i].disFn_fpts,eles[i].dA_fpts);
    }
  }
  else {
    /* Extrapolate transformed normal flux */
#pragma omp parallel for
    for (uint i=0; i<eles.size(); i++) {
      opers[eles[i].eType][eles[i].order].applyExtrapolateFn(eles[i].F_spts,eles[i].tNorm_fpts,eles[i].disFn_fpts);
    }
  }
}

void solver::correctDivFlux(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].calcDeltaFn();
    opers[eles[i].eType][eles[i].order].applyCorrectDivF(eles[i].dFn_fpts,eles[i].divF_spts[step]);
  }

}

void solver::calcGradU_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyGradSpts(eles[i].U_spts,eles[i].dU_spts);
  }
}

void solver::correctU()
{

}

void solver::extrapolateGradU()
{

}

void solver::calcEntropyErr_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].calcEntropyErr_spts();
  }
}

double solver::calcNormError(int field)
{
  double err = 0;
#pragma omp parallel for reduction(+:err)
  for (uint i=0; i<eles.size(); i++) {
    auto tmpErr = eles[i].getNormError();
    err += tmpErr[field];
  }

  return sqrt(err);
}

void solver::moveMesh(int step)
{
  if (!params->motion) return;

#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].move(step);
  }
}

void solver::setupOperators()
{
  // Get all element types & olynomial orders in mesh
  for (auto& e:eles) {
    eTypes.insert(e.eType);
    polyOrders[e.eType].insert(e.order);
  }

  for (auto& e: eTypes) {
    for (auto& p: polyOrders[e]) {
      opers[e][p].setupOperators(e,p,Geo,params);
    }
  }
}

void solver::setupElesFaces(void) {
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].setup(params,Geo);
  }

  // Finish setting up internal & boundary faces
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->setupFace();
  }

  // Finish setting up MPI faces
#pragma omp parallel for
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->setupFace();
  }
}

void solver::finishMpiSetup(void)
{
#pragma omp parallel for
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->finishRightSetup();
  }
}

void solver::readRestartFile(void) {

  ifstream dataFile;
  dataFile.precision(15);

  // Get the file name & open the file
  char fileNameC[50];
  string fileName = params->dataFileName;
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],params->restartIter);

  dataFile.open(fileNameC);

  if (!dataFile.is_open())
    FatalError("Cannont open restart file.");

  // Find the start of the UnstructuredData region
  bool found = false;
  string str;
  while (getline(dataFile,str)) {
    stringstream ss;
    ss.str(str);
    ss >> str;
    if (str.compare("<UnstructuredGrid>")==0) {
      found = true;
      break;
    }
  }

  if (!found)
    FatalError("Cannot find UnstructuredData tag in restart file.");

  // Read restart data & setup all data arrays
  for (auto& e:eles) {
    e.restart(dataFile,params,Geo);
  }

  dataFile.close();

  // Setup all transformations and other geometry-related arrays
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].setupAllGeometry();
  }

  // Finish setting up internal faces
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->setupFace();
  }
}


void solver::readReferenceSolution(string fileName, int iter) {

  ifstream dataFile;
  dataFile.precision(15);

  // Get the file name & open the file
  char fileNameC[50];
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],iter);

  dataFile.open(fileNameC);

  if (!dataFile.is_open())
    FatalError("Cannont open reference-solution file.");

  // Find the start of the UnstructuredData region
  bool found = false;
  string str;
  while (getline(dataFile,str)) {
    stringstream ss;
    ss.str(str);
    ss >> str;
    if (str.compare("<UnstructuredGrid>")==0) {
      found = true;
      break;
    }
  }

  if (!found)
    FatalError("Cannot find UnstructuredData tag in reference-data file.");

  // Read restart data & setup all data arrays
  for (auto& e:eles) {
    e.readReferenceSolution(dataFile,params,Geo);
  }

  dataFile.close();
}

void solver::initializeSolution()
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].setInitialCondition();
  }
}

solver::~solver()
{
  // Clean up memory allocated with 'new'
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->~face();
  }

  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->~face();
  }
}
