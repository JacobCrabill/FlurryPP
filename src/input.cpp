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

void fileReader::openFile(string fileName) {
  this->fileName = fileName;
  optFile.open(fileName.c_str(), ifstream::in);
}

void fileReader::closeFile()
{
  optFile.close();
}

template<typename T>
void fileReader::getScalarValue(string optName, T &opt, T defaultVal)
{
  if (!optFile.is_open()) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given string
  size_t start, end;
  string str;
  while (!optFile.eof()) {
    getline(optFile,str);

    // Remove any leading whitespace & search for option string
    start = str.find_first_not_of(" ");
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

      return;
    }
  }

  opt = defaultVal;
}

template<typename T>
void fileReader::getScalarValue(string optName, T &opt)
{
  if (!optFile.is_open()) {
    optFile.open(fileName.c_str());
    if (!optFile.is_open())
      FatalError("Cannont open input file for reading.");
  }

  // Rewind to the start of the file
  optFile.seekg(0,optFile.beg);

  // Search for the given string
  size_t start, end;
  string str;
  while (!optFile.eof()) {
    getline(optFile,str);

    // Remove any leading whitespace & search for option string
    start = str.find_first_not_of(" ");
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

input::input(const input& inInput)
{

}

void input::readInputFile(char *filename)
{
  /* --- Open Input File --- */
  string fName;
  fName.assign(filename);

  opts.openFile(fName);

  /* --- Read input file & store all simulation parameters --- */

  opts.getScalarValue("equation",equation,1);
  opts.getScalarValue("viscous",viscous,0);
  opts.getScalarValue("order",order,3);
  opts.getScalarValue("riemann_type",riemann_type,0);
  opts.getScalarValue("ic_type",ic_type,0);
  opts.getScalarValue("test_case",test_case,0);

  opts.getScalarValue("iterMax",iterMax);
  opts.getScalarValue("restart",restart,0);
  opts.getScalarValue("restartIter",restartIter);

  opts.getScalarValue("mesh_type",mesh_type,0);
  opts.getScalarValue("mesh_file_name",meshFileName);
  opts.getScalarValue("nx",nx,10);
  opts.getScalarValue("ny",ny,10);
  opts.getScalarValue("xmin",xmin,-10.);
  opts.getScalarValue("xmax",xmax,10.);
  opts.getScalarValue("ymin",ymin,-10.);
  opts.getScalarValue("ymax",ymax,10.);

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
}
