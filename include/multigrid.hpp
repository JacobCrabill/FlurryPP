/*!
 * \file multigrid.hpp
 * \brief Class to apply h/p multigrid operations
 *
 * \author - Josh Romero
 *           Adapted for Flurry++ by Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2016 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include <memory>
#include <vector>

#include "input.hpp"
#include "matrix.hpp"
#include "solver.hpp"


class multiGrid
{
  private:
    input *params = NULL;
    vector<input> pInputs, hInputs;
    int order;
    vector<shared_ptr<solver>> pGrids, hGrids;
    vector<shared_ptr<geo>> pGeos, hGeos;
    shared_ptr<geo> fine_grid;

    //! For nested HMG method
    vector<vector<int>> parent_cells;
    vector<matrix<int>> child_cells;

    void restrict_pmg(solver &grid_fine, solver &grid_coarse);
    void prolong_err(solver &grid_c, solver &grid_f);
    void compute_source_term(solver &grid);

    void restrict_hmg(solver &grid_f, solver&grid_c, uint H);
    void prolong_hmg(solver &grid_c, solver&grid_f, uint H);
    void setup_h_level(geo& mesh_c, geo& mesh_f, int H, int refine_level);

  public:
    void setup(int order, input *params, solver& Solver);
    void cycle(solver &Solver);
};
