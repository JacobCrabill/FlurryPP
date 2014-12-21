/*!
 * \file geo.cpp
 * \brief Class for handling geometry setup & modification
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

#include "../include/geo.hpp"


geo::setup(input* _params)
{
  params = _params;

  switch(params->mesh_type) {
  case (READ_MESH):
    readGmsh(params->mesh_file_name);

  case (CREATE_MESH):

  default:
    FatalError("Mesh type not recognized.")
  }
}
