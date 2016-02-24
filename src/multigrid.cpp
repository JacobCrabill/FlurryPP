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
#include "output.hpp"
#include "solver.hpp"

void multiGrid::setup(int order, input *params, solver &Solver)
{
  this->order = order;
  this->params = params;

  pInputs.assign(order, *params);
  pGrids.resize(order);

  /* H-P Multigrid using Refinement Method */
  if (params->HMG)
  {
    if (params->lowOrder != 0)
      FatalError("H-Multigrid only supported for PMG lowOrder = 0.");

    if (params->nDims == 3)
      FatalError("H-Multigrid only supported for 2D currently.");

    geo coarse_grid;
    coarse_grid.setup(params,true);

    hInputs.assign(params->n_h_levels, *params);
    hGrids.resize(params->n_h_levels);
    hGeos.resize(params->n_h_levels);
    pGeos.resize(params->order);

    for (int H = 0; H < params->n_h_levels; H++)
    {
      if (params->rank == 0) cout << endl << "H-Multigrid: Setting up H = " << H << endl;

      hInputs[H].dataFileName += "_H" + std::to_string(H) + "_";
      hGrids[H] = make_shared<solver>();
      hGeos[H] = make_shared<geo>();

      /* Refine the initial coarse grid to produce the fine grids */
      setup_h_level(coarse_grid, *hGeos[H], params->n_h_levels - H - 1);

      hGrids[H]->setup(&hInputs[H], params->lowOrder, &(*hGeos[H]));
      hGrids[H]->initializeSolution(true);
    }

    /* Create final fine grid and re-setup the given solver */
    if (params->rank == 0) cout << endl << "H-Multigrid: Setting up fine grid solver" << endl;
    fine_grid = make_shared<geo>();
    setup_h_level(coarse_grid, *fine_grid, params->n_h_levels);
    Solver.setup(params, order, &(*fine_grid));
    Solver.initializeSolution();

    /* Instantiate P-grid solvers using finest mesh */
    for (int P = 0; P < order; P++)
    {
      if (P < params->lowOrder)
      {
        pGrids[P] = NULL;
      }
      else
      {
        if (params->rank == 0) cout << endl << "P-Multigrid: Setting up P = " << P << endl;

        pInputs[P].dataFileName += "_P" + std::to_string(P) + "_";
        pGeos[P] = make_shared<geo>(*fine_grid);
        pGrids[P] = make_shared<solver>();

        pGeos[P]->setup_hmg(&pInputs[P], fine_grid->gridID, fine_grid->gridRank, fine_grid->nProcGrid, fine_grid->gridIdList);
        pGrids[P]->setup(&pInputs[P], P, &(*pGeos[P]));
        pGrids[P]->initializeSolution(true);
      }
    }
  }

  /* P-Multigrid Alone */
  else {
    /* Instantiate coarse grid solvers */
    for (int P = 0; P < order; P++)
    {
      if (P < params->lowOrder)
      {
        pGrids[P] = NULL;
      }
      else
      {
        if (params->rank == 0) cout << endl << "P-Multigrid: Setting up P = " << P << endl;

        pInputs[P].dataFileName += "_P" + std::to_string(P) + "_";
        pGrids[P] = make_shared<solver>();
        pGrids[P]->setup(&pInputs[P], P);
        pGrids[P]->initializeSolution(true);
      }
    }

    if (params->rank == 0) cout << endl << "P-Multigrid: Setting up P = " << params->order << endl;
    Solver.setup(params, params->order);
    Solver.initializeSolution();
  }

  /* Still some weird bug in initialization of PMG solvers; this is a
   * workaround for the moment */
  cycle(Solver);
  Solver.initializeSolution();
}

void multiGrid::setup_h_level(geo &mesh_c, geo &mesh_f, int refine_level)
{
  if (refine_level > 0)
    refineGrid2D(mesh_c, mesh_f, refine_level, mesh_c.nNodesPerCell, params->shapeOrder);
  else
    mesh_f = mesh_c;

  vector<int> epart(0);
#ifndef _NO_MPI
  if (mesh_c.nproc > 1) {
    int nSplit = 1 << params->nDims;
    nSplit = std::pow(nSplit, refine_level);

    epart.resize(mesh_f.nEles);
    for (uint e = 0; e < mesh_f.nEles; e++) {
      epart[e] = mesh_c.epart[e/nSplit];
    }
  }
#endif

  mesh_f.setup_hmg(params, mesh_c.gridID, mesh_c.gridRank, mesh_c.nProcGrid, mesh_c.gridIdList, epart);
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

    if (P-1 >= (int) params->lowOrder || params->HMG)
    {
      /* Update residual and add source */
      pGrids[P]->calcResidual(0);

#pragma omp parallel for
      for (uint e = 0; e < pGrids[P]->eles.size(); e++)
      {
        pGrids[P]->eles[e]->divF_spts[0] += pGrids[P]->eles[e]->src_spts;
      }

      if (P-1 >= (int) params->lowOrder)
      {
        /* Restrict to next coarse grid */
        restrict_pmg(*pGrids[P], *pGrids[P-1]);
      }
    }
  }

  if (params->HMG)
  {
    /* Initial restriction to first HMG level */
    restrict_hmg(*pGrids[params->lowOrder], *hGrids[0], 0);

    for (int H = 0; H < params->n_h_levels; H++)
    {
      /* Generate source term */
      compute_source_term(*hGrids[H]);

      /* Copy initial solution to solution storage */
  #pragma omp parallel for
      for (uint e = 0; e < hGrids[H]->eles.size(); e++)
      {
        hGrids[H]->eles[e]->sol_spts = hGrids[H]->eles[e]->U_spts;
      }

      /* Update solution on coarse level */
      for (uint step = 0; step < params->smoothSteps; step++)
      {
        hGrids[H]->update(true);
      }

      if (H+1 < params->n_h_levels)
      {
        /* Update residual and add source */
        hGrids[H]->calcResidual(0);

#pragma omp parallel for
        for (uint e = 0; e < hGrids[H]->eles.size(); e++)
        {
          hGrids[H]->eles[e]->divF_spts[0] += hGrids[H]->eles[e]->src_spts;
        }

        /* Restrict to next coarse grid */
        restrict_hmg(*hGrids[H], *hGrids[H+1], H+1);
      }
    }

    /* --- Upward HMG Cycle --- */
    for (int H = params->n_h_levels-1; H >= 0; H--)
    {
      /* Advance again (v-cycle)*/
      for (unsigned int step = 0; step < params->smoothSteps; step++)
      {
        hGrids[H]->update(true);
      }

      /* Generate error */
  #pragma omp parallel for
      for (int e = 0; e < hGrids[H]->eles.size(); e++)
      {
        hGrids[H]->eles[e]->corr_spts  = hGrids[H]->eles[e]->U_spts;
        hGrids[H]->eles[e]->corr_spts -= hGrids[H]->eles[e]->sol_spts;
      }

      /* Prolong error and add to fine grid solution */
      if (H > 0)
      {
        prolong_hmg(*hGrids[H], *hGrids[H-1], H-1);
      }
      else
      {
        prolong_hmg(*hGrids[0], *pGrids[params->lowOrder], 0);
      }
    }
  }

  /* --- Upward PMG Cycle --- */
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

void multiGrid::restrict_hmg(solver &grid_f, solver &grid_c, uint H)
{
  int nSplit = 1 << params->nDims;

#pragma omp parallel for
  for (uint e = 0; e < grid_c.eles.size(); e++)
  {
    auto &e_c = *grid_c.eles[e];
    e_c.U_spts.initializeToZero();
    e_c.divF_spts[0].initializeToZero();

    for (uint j = 0; j < nSplit; j++)
    {
      // what to do...? Take avg (Josh's simple method), or do Galerkin projection?
      // For the moment: take avg value [Only remotely reasonable for P = 0]
      auto &e_f = *grid_f.eles[e*nSplit+j];
      e_c.U_spts += e_f.U_spts;
      e_c.divF_spts[0] += e_f.divF_spts[0];
    }
    e_c.U_spts /= nSplit;
    e_c.divF_spts[0] /= nSplit;
  }
}

void multiGrid::prolong_hmg(solver &grid_c, solver &grid_f, uint H)
{
  int nSplit = 1 << params->nDims;

#pragma omp parallel for
  for (uint e = 0; e < grid_f.eles.size(); e++)
  {
    auto &e_f = *grid_f.eles[e];
    auto &e_c = *grid_c.eles[e/nSplit];
    // what to do...? Copy value (Josh's simple method), or do Galerkin projection?
    // For the moment: simple method [Only remotely reasonable for P = 0]
    e_f.U_spts += e_c.corr_spts;
  }
}
