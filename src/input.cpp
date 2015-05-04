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

#include <fstream>
#include <sstream>
#include <string>

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
  FatalError(errMsg.c_str())
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
    FatalError(errMsg.c_str())
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
  FatalError(errMsg.c_str())
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
  opts.getScalarValue("ic_type",ic_type,0);
  if (equation==ADVECTION_DIFFUSION) {
    opts.getScalarValue("advectVx",advectVx,1.);
    opts.getScalarValue("advectVy",advectVy,1.);
    opts.getScalarValue("lambda",lambda,1.);
    nFields = 1;
  } else if (equation==NAVIER_STOKES) {
    if (ic_type == 0) {
      opts.getScalarValue("rhoIC",rhoIC,1.);
      opts.getScalarValue("vxIC",vxIC,.2);
      opts.getScalarValue("vyIC",vyIC,0.);
      opts.getScalarValue("pIC",pIC,.7142857143);
    }
    opts.getScalarValue("gamma",gamma,1.4);
    opts.getScalarValue("RGas",RGas,286.9);
    opts.getScalarValue("rhoBound",rhoBound,1.);
    opts.getScalarValue("uBound",uBound,.2);
    opts.getScalarValue("vBound",vBound,0.);
    opts.getScalarValue("wBound",wBound,0.);
    opts.getScalarValue("pBound",pBound,.7142857143);
    opts.getScalarValue("TBound",TBound,300.);
    opts.getScalarValue("TWall",TWall,300.);
    opts.getScalarValue("beta",beta,2.);
    opts.getScalarValue("slipPenalty",slipPenalty,false);
    nFields = 4;
  }

  opts.getScalarValue("dt",dt);
  opts.getScalarValue("viscous",viscous,0);
  opts.getScalarValue("motion",motion,0);
  opts.getScalarValue("order",order,3);
  opts.getScalarValue("riemann_type",riemann_type,0);
  opts.getScalarValue("test_case",test_case,0);
  opts.getScalarValue("iterMax",iterMax);

  opts.getScalarValue("timeType",timeType,0);

  opts.getScalarValue("restart",restart,0);
  if (restart) {
    opts.getScalarValue("restartIter",restartIter);
  }

  opts.getScalarValue("mesh_type",mesh_type,1); // CREATE_MESH by default

  if (mesh_type == CREATE_MESH) {
    opts.getScalarValue("nDims",nDims,2);
    opts.getScalarValue("nx",nx,10);
    opts.getScalarValue("ny",ny,10);
    opts.getScalarValue("xmin",xmin,-10.);
    opts.getScalarValue("xmax",xmax,10.);
    opts.getScalarValue("ymin",ymin,-10.);
    opts.getScalarValue("ymax",ymax,10.);
    opts.getScalarValue("create_bcTop",create_bcTop,string("periodic"));
    opts.getScalarValue("create_bcBottom",create_bcBottom,string("periodic"));
    opts.getScalarValue("create_bcLeft",create_bcLeft,string("periodic"));
    opts.getScalarValue("create_bcRight",create_bcRight,string("periodic"));
  }else if (mesh_type == READ_MESH) {
    opts.getScalarValue("mesh_file_name",meshFileName);
    opts.getMap("mesh_bound",meshBounds);
  }

  opts.getScalarValue("monitor_res_freq",monitor_res_freq,10);
  opts.getScalarValue("resType",resType,2);
  opts.getScalarValue("plot_freq",plot_freq,100);
  opts.getScalarValue("plot_type",plot_type,1);
  opts.getScalarValue("restart_freq",restart_freq,100);
  opts.getScalarValue("dataFileName",dataFileName,string("simData"));

  opts.getScalarValue("spts_type_tri",sptsTypeTri,string("Legendre"));
  opts.getScalarValue("spts_type_quad",sptsTypeQuad,string("Legendre"));
  opts.getScalarValue("vcjhSchemeTri",vcjhSchemeTri,0);
  opts.getScalarValue("vcjhSchemeQuad",vcjhSchemeQuad,0);

  /* --- Cleanup ---- */
  opts.closeFile();

  /* --- Additional Processing --- */
  if (restart) {
    initIter = restartIter;
  }else{
    initIter = 0;
  }

  iter = initIter;
}
