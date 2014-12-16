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

#include <map>

#include "global.hpp"
#include "ele.hpp"
#include "operators.hpp"

class solver
{
public:
  //! Map from eType to element-specific operator
  map<int,oper> opers;

  //!
  vector<ele> eles;

private:
}

