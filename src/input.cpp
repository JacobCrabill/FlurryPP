/*!
 * \file input.cpp
 * \brief Class to read & store simulation parameters
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
#include <sstream>

fileReader::fileReader()
{

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
  size_t start, end, diff;
  string str;

  openFile();

  if (!optFile.is_open() || !getline(optFile,str)) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given string
  while (getline(optFile,str)) {

    // Remove any leading whitespace & search for option string
    start = str.find_first_not_of(" ");
    diff = str.length() - start;
    if (diff<optName.length()) continue; // If line is shorter that length of opt string, skip

    if (!str.compare(start,optName.length(),optName)) {

      // Find start/end of scalar value
      start += optName.length()+1;
      end = str.find_first_of(" ");
      if (end==string::npos) {
        end = str.length();
      }

      // Put scalar into stringstream for conversion
      if (!(stringstream(str.substr(start,end-start)) >> opt)) {
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
  size_t start, end, diff;
  string str;

  openFile();

  if (!optFile.is_open()) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given string  
  while (getline(optFile,str)) {

    // Remove any leading whitespace & search for option string
    start = str.find_first_not_of(" ");
    diff = str.length() - start;
    if (diff<optName.length()) continue; // If line is shorter that length of opt string, skip

    if (!str.compare(start,optName.length(),optName)) {

      // Find start/end of scalar value
      start += optName.length()+1;
      end = str.find_first_of(" ");
      if (end==string::npos) {
        end = str.length();
      }

      // Put scalar into stringstream for conversion
      if (!(stringstream(str.substr(start,end-start)) >> opt)) {
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
  if (equation==ADVECTION_DIFFUSION) {
    opts.getScalarValue("advectVx",advectVx,1.);
    opts.getScalarValue("advectVy",advectVy,1.);
  } else if (equation==NAVIER_STOKES) {

  }

  opts.getScalarValue("viscous",viscous,0);
  opts.getScalarValue("motion",motion,0);
  opts.getScalarValue("order",order,3);
  opts.getScalarValue("riemann_type",riemann_type,0);
  opts.getScalarValue("ic_type",ic_type,0);
  opts.getScalarValue("test_case",test_case,0);
  opts.getScalarValue("iterMax",iterMax);

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
  }else if (mesh_type == READ_MESH) {
    opts.getScalarValue("mesh_file_name",meshFileName);
  }

  opts.getScalarValue("plot_freq",plot_freq,100);
  opts.getScalarValue("restart_freq",restart_freq,100);
  opts.getScalarValue("dataFileName",dataFileName,string("simData"));

  opts.getScalarValue("spts_type_tri",sptsTypeTri,string("Legendre"));
  opts.getScalarValue("spts_type_quad",sptsTypeQuad,string("Legendre"));

  /* --- Cleanup ---- */
  opts.closeFile();

  /* --- Additional Processing --- */
  if (restart) {
    initIter = restartIter;
  }else{
    initIter = 0;
  }

  // --- Testing ---
  cout << "equation = " << equation << endl;
  cout << "viscous = " << viscous << endl;
  cout << "order = " << order << endl;
  cout << "spts_type_tri = " << sptsTypeTri << endl;
  cout << "mesh_type = " << mesh_type<< endl;

}
