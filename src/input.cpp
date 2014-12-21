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


void fileReader::getScalarOpt(string optName, fileReader::T &opt, fileReader::T defaultVal)
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

  opts.getScalarValue("mesh_type",mesh_type,0);
  opts.getScalarValue("mesh_file_name",mesh_file_name);
  opts.getScalarValue("nx",nx,10);
  opts.getScalarValue("ny",ny,10);
  opts.getScalarValue("xmin",xmin,-10);
  opts.getScalarValue("xmax",xmax,10);
  opts.getScalarValue("ymin",ymin,-10);
  opts.getScalarValue("ymax",ymax,10);

  opts.getScalarValue("spts_type_tri",spts_type_tri,"Legendre");
  opts.getScalarValue("spts_type_quad",spts_type_quad,"Legendre");


  // Basic value-getter:
  //opts.getScalarValue("",);

  /* --- Cleanup ---- */
  opts.closeFile();
}
