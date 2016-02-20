/*!
 * \file multigrid.cpp
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

#include "multigrid.hpp"

#include <sstream>
#include <omp.h>

#include "input.hpp"
#include "solver.hpp"

void multiGrid::setup(int order, input *params, solver &Solver)
{
  this->order = order;
  this->params = params;

  Inputs.resize(order);

  if (params->HMG)
  {
    Inputs.push_back(*params);
    Inputs.back().order = 0;
    for (int H = 0; H < params->n_h_levels; H++)
    {
      hGrids.push_back(make_shared<solver>());
      hGeos.push_back(make_shared<geo>());

      /* Refine the initial coarse grid to produce the fine grids */
      setup_h_level(*Solver.Geo, *hGeos.back(), H+1);

      hGrids[H]->setup(&Inputs.back(), &(*hGeos.back()));

    }

    /* Swap the initial (coarse-grid) solver for the new fine-grid solver */
    pGrids.push_back(make_shared<solver>(Solver)); // does this do what I want it to do?
    Solver = *hGrids.back();
  }

  /* Instantiate coarse grid solvers */
  for (int P = 0; P < order; P++)
  {
    if (P < params->lowOrder)
    {
      pGrids.push_back(NULL);
    }
    else
    {
      if (params->rank == 0) std::cout << "P = " << P << std::endl;
      Inputs[P] = *params;
      Inputs[P].order = P;
      pGrids.push_back(make_shared<solver>());
      pGrids[P]->setup(&Inputs[P]);
      pGrids[P]->initializeSolution(true);
    }
  }
}

void multiGrid::setup_h_level(geo &mesh_c, geo &mesh_f, int level)
{
  /* use shape_quad, mesh_c.xv/c2v/e2c?/etc. to split all edges, add eles &
   * points, etc.
   * Add param 'shape_order' to control order of fine-grid elements [assuming
   * coarse grid is of high-order] */
}

void multiGrid::cycle(solver &Solver)
{
  /* Update residual on finest grid level and restrict */
  Solver.calcResidual(0);

  restrict_pmg(Solver, *pGrids[order-1]);

  for (int P = order-1; P >= (int) params->lowOrder; P--)
  {
    /* Generate source term */
    compute_source_term(*pGrids[P]);

    /* Copy initial solution to solution storage */
#pragma omp parallel for
    for (uint e = 0; e < pGrids[P]->eles.size(); e++)
    {
      pGrids[P]->eles[e]->sol_spts = pGrids[P]->eles[e]->U_spts;
    }

    /* Update solution on coarse level */
    for (uint step = 0; step < params->smoothSteps; step++)
    {
      pGrids[P]->update(true);
    }

    if (P-1 >= (int) params->lowOrder)
    {
      /* Update residual and add source */
      pGrids[P]->calcResidual(0);

#pragma omp parallel for
      for (uint e = 0; e < pGrids[P]->eles.size(); e++)
      {
        pGrids[P]->eles[e]->divF_spts[0] += pGrids[P]->eles[e]->src_spts;
      }

      /* Restrict to next coarse grid */
      restrict_pmg(*pGrids[P], *pGrids[P-1]);
    }
  }

  for (int P = (int) params->lowOrder; P <= order-1; P++)
  {
    /* Advance again (v-cycle)*/
    for (unsigned int step = 0; step < params->smoothSteps; step++)
    {
      pGrids[P]->update(true);
    }

    /* Generate error */
#pragma omp parallel for
    for (int e = 0; e < pGrids[P]->eles.size(); e++)
    {
      pGrids[P]->eles[e]->corr_spts  = pGrids[P]->eles[e]->U_spts;
      pGrids[P]->eles[e]->corr_spts -= pGrids[P]->eles[e]->sol_spts;
    }

    /* Prolong error and add to fine grid solution */
    if (P < order-1)
    {
      prolong_err(*pGrids[P], *pGrids[P+1]);
    }
  }

  /* Prolong correction and add to finest grid solution */
  prolong_err(*pGrids[order-1], Solver);
}

void multiGrid::restrict_pmg(solver &grid_f, solver &grid_c)
{
  if (grid_f.order - grid_c.order > 1)
    FatalError("Cannot restrict more than 1 order currently!");

#pragma omp parallel for
  for (uint e = 0; e < grid_f.eles.size(); e++)
  {
    auto &e_f = *grid_f.eles[e];
    auto &e_c = *grid_c.eles[e];
    auto &opp_res = grid_f.opers[e_f.eType][e_f.order].opp_restrict;

    /* Restrict solution */
    opp_res.timesMatrix(e_f.U_spts, e_c.U_spts);

    /* Restrict residual */
    opp_res.timesMatrix(e_f.divF_spts[0], e_c.divF_spts[0]);
  }
}

void multiGrid::prolong_pmg(solver &grid_c, solver &grid_f)
{
  if (grid_f.order - grid_c.order > 1)
    FatalError("Cannot prolong more than 1 order currently!");

#pragma omp parallel for
  for (uint e = 0; e < grid_c.eles.size(); e++)
  {
    auto &e_c = *grid_c.eles[e];
    auto &e_f = *grid_f.eles[e];
    auto &opp_pro = grid_c.opers[e_c.eType][e_c.order].opp_prolong;
    opp_pro.timesMatrix(e_c.U_spts, e_f.U_spts);
  }
}

void multiGrid::prolong_err(solver &grid_c, solver &grid_f)
{
#pragma omp parallel for
  for (uint e = 0; e < grid_c.eles.size(); e++)
  {
    auto &e_c = *grid_c.eles[e];
    auto &e_f = *grid_f.eles[e];
    auto &opp_pro = grid_c.opers[e_c.eType][e_c.order].opp_prolong;
    opp_pro.timesMatrixPlus(e_c.corr_spts, e_f.U_spts);
  }
}


void multiGrid::compute_source_term(solver &grid)
{
  /* Copy restricted fine grid residual to source term */
#pragma omp parallel for
  for (uint e = 0; e < grid.eles.size(); e++)
  {
    grid.eles[e]->src_spts = grid.eles[e]->divF_spts[0];
  }

  /* Update residual on current coarse grid */
  grid.calcResidual(0);

  /* Subtract to generate source term */
#pragma omp parallel for
  for (uint e = 0; e < grid.eles.size(); e++)
  {
    grid.eles[e]->src_spts -= grid.eles[e]->divF_spts[0];
  }
}

void multiGrid::restrict_hmg(solver &grid_f, solver&grid_c)
{

}

void multiGrid::prolong_hmg(solver &grid_c, solver&grid_f)
{

}
