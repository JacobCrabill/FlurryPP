/*!
 * \file input.cpp
 * \brief Class to read & store simulation parameters from file
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Fux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include "../include/input.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>

fileReader::fileReader()
{

}

fileReader::fileReader(string fileName)
{
  this->fileName = fileName;
}

fileReader::~fileReader()
{
  if (optFile.is_open()) optFile.close();
}

void fileReader::setFile(string fileName)
{
  this->fileName = fileName;
}

void fileReader::openFile(void)
{
  optFile.open(fileName.c_str(), ifstream::in);
}

void fileReader::closeFile()
{
  optFile.close();
}

template<typename T>
void fileReader::getScalarValue(string optName, T &opt, T defaultVal)
{
  string str, optKey;

  openFile();

  if (!optFile.is_open() || !getline(optFile,str)) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given option string
  while (getline(optFile,str)) {
    // Remove any leading whitespace & see if first word is the input option
    stringstream ss;
    ss.str(str);
    ss >> optKey;
    if (optKey.compare(optName)==0) {
      if (!(ss >> opt)) {
        // This could happen if, for example, trying to assign a string to a double
        cout << "WARNING: Unable to assign value to option " << optName << endl;
        cout << "Using default value of " << defaultVal << " instead." << endl;
        opt = defaultVal;
      }

      closeFile();
      return;
    }
  }

  opt = defaultVal;
  closeFile();
}

template<typename T>
void fileReader::getScalarValue(string optName, T &opt)
{
  string str, optKey;

  openFile();

  if (!optFile.is_open()) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given option string
  while (getline(optFile,str)) {
    // Remove any leading whitespace & see if first word is the input option
    stringstream ss;
    ss.str(str);
    ss >> optKey;
    if (optKey.compare(optName)==0) {
      if (!(ss >> opt)) {
        // This could happen if, for example, trying to assign a string to a double
        cerr << "WARNING: Unable to assign value to option " << optName << endl;
        string errMsg = "Required option not set: " + optName;
        FatalError(errMsg.c_str())
      }

      closeFile();
      return;
    }
  }

  // Option was not found; throw error & exit
  string errMsg = "Required option not found: " + optName;
  FatalError(errMsg.c_str());
}

template<typename T, typename U>
void fileReader::getMap(string optName, map<T,U> &opt) {
  string str, optKey;
  T tmpT;
  U tmpU;
  bool found;

  openFile();

  if (!optFile.is_open()) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given option string
  while (getline(optFile,str)) {
    // Remove any leading whitespace & see if first word is the input option
    stringstream ss;
    ss.str(str);
    ss >> optKey;
    if (optKey.compare(optName)==0) {
      found = true;
      if (!(ss >> tmpT >> tmpU)) {
        // This could happen if, for example, trying to assign a string to a double
        cerr << "WARNING: Unable to assign value to option " << optName << endl;
        string errMsg = "Required option not set: " + optName;
        FatalError(errMsg.c_str())
      }

      opt[tmpT] = tmpU;
      optKey = "";
    }
  }

  if (!found) {
    // Option was not found; throw error & exit
    string errMsg = "Required option not found: " + optName;
    FatalError(errMsg.c_str());
  }

  closeFile();
}

template<typename T>
void fileReader::getVectorValue(string optName, vector<T> &opt)
{
  string str, optKey;

  openFile();

  if (!optFile.is_open()) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given option string
  while (getline(optFile,str)) {
    // Remove any leading whitespace & see if first word is the input option
    stringstream ss;
    ss.str(str);
    ss >> optKey;
    if (optKey.compare(optName)==0) {
      int nVals;
      if (!(ss >> nVals)) {
        // This could happen if, for example, trying to assign a string to a double
        cerr << "WARNING: Unable to read number of entries for vector option " << optName << endl;
        string errMsg = "Required option not set: " + optName;
        FatalError(errMsg.c_str());
      }

      opt.resize(nVals);
      for (int i=0; i<nVals; i++) {
        if (!ss >> opt[i]) {
          cerr << "WARNING: Unable to assign all values to vector option " << optName << endl;
          string errMsg = "Required option not set: " + optName;
          FatalError(errMsg.c_str())
        }
      }

      closeFile();
      return;
    }
  }

  // Option was not found; throw error & exit
  string errMsg = "Required option not found: " + optName;
  FatalError(errMsg.c_str());
}

input::input()
{

}

void input::readInputFile(char *filename)
{
  /* --- Open Input File --- */
  string fName;
  fName.assign(filename);

  opts.setFile(fName);

  /* --- Read input file & store all simulation parameters --- */

  opts.getScalarValue("equation",equation,1);
  opts.getScalarValue("icType",icType,0);
  if (equation==ADVECTION_DIFFUSION) {
    opts.getScalarValue("advectVx",advectVx,1.);
    opts.getScalarValue("advectVy",advectVy,1.);
    opts.getScalarValue("advectVz",advectVz,0.);
    opts.getScalarValue("diffD",diffD,1.);
    opts.getScalarValue("lambda",lambda,1.);
    nFields = 1;
  }
  else if (equation==NAVIER_STOKES) {
    opts.getScalarValue("gamma",gamma,1.4);
    opts.getScalarValue("RGas",RGas,286.9);
    opts.getScalarValue("rhoBound",rhoBound,1.);
    opts.getScalarValue("uBound",uBound,.2);
    opts.getScalarValue("vBound",vBound,0.);
    opts.getScalarValue("wBound",wBound,0.);
    opts.getScalarValue("pBound",pBound,.7142857143);
    opts.getScalarValue("TBound",TBound,300.);
    opts.getScalarValue("TWall",TWall,300.);
    opts.getScalarValue("entropySensor",calcEntropySensor,false);
    opts.getScalarValue("slipPenalty",slipPenalty,false);
    if (slipPenalty) {
     opts.getScalarValue("Kp",Kp);
     opts.getScalarValue("Kd",Kd);
     opts.getScalarValue("Ki",Ki);
    }
    if (icType == 0) {
      opts.getScalarValue("rhoIC",rhoIC,rhoBound);
      opts.getScalarValue("vxIC",vxIC,uBound);
      opts.getScalarValue("vyIC",vyIC,vBound);
      opts.getScalarValue("vzIC",vzIC,wBound);
      opts.getScalarValue("pIC",pIC,pBound);
    }
  }

  opts.getScalarValue("timeType",timeType,4);
  opts.getScalarValue("dtType",dtType,0);
  if (dtType == 1)
    opts.getScalarValue("CFL",CFL);
  else
    opts.getScalarValue("dt",dt);

  opts.getScalarValue("viscous",viscous,0);
  opts.getScalarValue("motion",motion,0);
  opts.getScalarValue("order",order,3);
  opts.getScalarValue("riemannType",riemannType,0);
  opts.getScalarValue("testCase",test_case,0);
  opts.getScalarValue("iterMax",iterMax);

  if (viscous && equation == NAVIER_STOKES) {
    opts.getScalarValue("Re",Re);
    opts.getScalarValue("Lref",Lref,1.0);
    opts.getScalarValue("muGas",muGas,1.827e-5);
    opts.getScalarValue("prandtl",prandtl,.72);
    opts.getScalarValue("SGas",SGas,120.);
    opts.getScalarValue("TGas",TGas,291.15);

    opts.getScalarValue("MachBound",MachBound,1.);
    opts.getScalarValue("nxBound",nxBound,1.);
    opts.getScalarValue("nyBound",nyBound,0.);
    opts.getScalarValue("nzBound",nzBound,0.);

    /* --- LDG Flux Parameters --- */
    opts.getScalarValue("LDG_penFact",penFact,0.5);
    opts.getScalarValue("LDG_tau",tau,1.0);
  }

  opts.getScalarValue("restart",restart,0);
  if (restart) {
    opts.getScalarValue("restartIter",restartIter);
  }

  opts.getScalarValue("meshType",meshType);

  if (meshType == CREATE_MESH) {
    opts.getScalarValue("nDims",nDims,2);
    opts.getScalarValue("nx",nx,10);
    opts.getScalarValue("ny",ny,10);
    opts.getScalarValue("nz",nz,10);
    opts.getScalarValue("xmin",xmin,-10.);
    opts.getScalarValue("xmax",xmax,10.);
    opts.getScalarValue("ymin",ymin,-10.);
    opts.getScalarValue("ymax",ymax,10.);
    opts.getScalarValue("zmin",zmin,-10.);
    opts.getScalarValue("zmax",zmax,10.);
    opts.getScalarValue("create_bcTop",create_bcTop,string("periodic"));
    opts.getScalarValue("create_bcBottom",create_bcBottom,string("periodic"));
    opts.getScalarValue("create_bcLeft",create_bcLeft,string("periodic"));
    opts.getScalarValue("create_bcRight",create_bcRight,string("periodic"));
    opts.getScalarValue("create_bcFront",create_bcFront,string("periodic"));
    opts.getScalarValue("create_bcBack",create_bcBack,string("periodic"));
  }
  else if (meshType == READ_MESH) {
    opts.getScalarValue("meshFileName",meshFileName);
    // Get mesh boundaries, boundary conditions & convert to lowercase
    map<string,string> meshBndTmp;
    opts.getMap("mesh_bound",meshBndTmp);
    for (auto& B:meshBndTmp) {
      string tmp1, tmp2;
      tmp1 = B.first; tmp2 = B.second;
      std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
      std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
      meshBounds[tmp1] = tmp2;
    }
  }
  opts.getScalarValue("periodicDX",periodicDX,(double)INFINITY);
  opts.getScalarValue("periodicDY",periodicDY,(double)INFINITY);
  opts.getScalarValue("periodicDZ",periodicDZ,(double)INFINITY);
  opts.getScalarValue("periodicTol",periodicTol,1e-6);

  opts.getScalarValue("monitorResFreq",monitorResFreq,10);
  opts.getScalarValue("resType",resType,2);
  opts.getScalarValue("plotFreq",plotFreq,100);
  opts.getScalarValue("plotType",plotType,1);
  opts.getScalarValue("restart_freq",restart_freq,100);
  opts.getScalarValue("dataFileName",dataFileName,string("simData"));

  opts.getScalarValue("spts_type_tri",sptsTypeTri,string("Legendre"));
  opts.getScalarValue("spts_type_quad",sptsTypeQuad,string("Legendre"));
  opts.getScalarValue("vcjhSchemeTri",vcjhSchemeTri,0);
  opts.getScalarValue("vcjhSchemeQuad",vcjhSchemeQuad,0);

  /* --- Shock Capturing --- */
  opts.getScalarValue("shockCapture",scFlag,0);
  if(scFlag == 1)
    opts.getScalarValue("threshold",threshold,1.0);

  opts.getScalarValue("squeeze",squeeze,0);

  /* --- Cleanup ---- */
  opts.closeFile();

  /* --- Additional Processing --- */
  if (restart) {
    initIter = restartIter;
  }else{
    initIter = 0;
  }

  if (squeeze) {
    // Entropy bound for polynomial squeezing
    exps0 = .9*pBound/(gamma*rhoBound);
  }

  iter = initIter;

  if (viscous) nonDimensionalize();
}


void input::nonDimensionalize(void)
{
  /* --- Calculate Reference / Freestream Non-Dimensionalized Values --- */

  // Normalize the boundary flow direction
  double nMag = sqrt(nxBound*nxBound+nyBound*nyBound+nzBound*nzBound);
  nxBound /= nMag;
  nyBound /= nMag;
  nzBound /= nMag;

  double Tref = TBound;

  // Calculate total velocity & individual components
  double UBound = MachBound * sqrt(gamma*RGas*TBound);
  uBound = UBound * nxBound;
  vBound = UBound * nyBound;
  wBound = UBound * nzBound;

  muBound = muGas * pow(TBound/TGas, 1.5) * ((TGas+SGas) / (TBound+SGas));

  rhoBound = muBound*Re/(UBound*Lref);
  pBound = rhoBound * RGas * TBound;

  double rhoRef = rhoBound;
  double pRef = rhoRef*UBound*UBound;
  double muRef = rhoRef*UBound*Lref;
  //double timeRef = Lref / UBound;
  //double RRef = (RGas*Tref) / (UBound*UBound);

  c_sth = SGas / TGas; // Sutherland's Law parameter
  mu_inf = muGas / muRef;
  rt_inf = TGas * RGas / (UBound*UBound);

  // Set up the dimensionless conditions at free-stream boundaries

  rhoBound = 1.;
  uBound = uBound / UBound;
  vBound = vBound / UBound;
  wBound = wBound / UBound;
  pBound = pBound / pRef;
  TtBound = (TBound/Tref) * (1. + 0.5*(gamma-1.)*MachBound*MachBound);
  PtBound = pBound * pow(1. + 0.5*(gamma-1.)*MachBound*MachBound, gamma/(gamma-1.));

  rhoIC = rhoBound;
  vxIC = uBound;
  vyIC = vBound;
  vzIC = wBound;
  pIC = pBound;
  muIC = muBound;
  TIC = TBound;
}
