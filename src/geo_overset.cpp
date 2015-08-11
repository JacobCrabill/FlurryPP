/*!
 * \file geo_overset.cpp
 * \brief Overset-related methods for the geo class
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill.
 *
 */

#include "geo.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>

#include "ele.hpp"
#include "face.hpp"
#include "overFace.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "mpi.h"
#include "metis.h"
#endif

void geo::registerGridDataTIOGA(void)
{
#ifndef _NO_MPI
  /* Note that this function should only be needed once during preprocessing */

  if (gridRank == 0)
    cout << "Geo: Grid " << gridID << ": Setting up TIOGA interface" << endl;

  // Allocate TIOGA grid processor
  tg = make_shared<tioga>();

  tg->setCommunicator(MPI_COMM_WORLD,params->rank,params->nproc);

  // Setup iwall, iover (nodes on wall & overset boundaries)
  iover.resize(0);
  iwall.resize(0);
  for (int ib=0; ib<nBounds; ib++) {
    if (bcList[ib] == OVERSET) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iover.push_back(bndPts(ib,iv));
      }
    }
    else if (bcList[ib] == SLIP_WALL || bcList[ib] == ADIABATIC_NOSLIP || bcList[ib] == ISOTHERMAL_NOSLIP) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iwall.push_back(bndPts(ib,iv));
      }
    }
  }

  int nwall = iwall.size();
  int nover = iover.size();
  int ntypes = 1;           //! Number of element types in grid block
  nodesPerCell = new int[1];
  nodesPerCell[0] = 8;      //! Number of nodes per element for each element type (but only one type so far)
  iblank.resize(nVerts);
  iblankCell.resize(nEles);
  iblankFace.resize(nFaces);

  if (c2v.getDim1() > 8) {
    // Quadratic elements present; setup a different c2v specifically for TIOGA
    tg_c2v.setup(nEles,8);
    for (int ic=0; ic<nEles; ic++) {
      for (int j=0; j<8; j++) {
        tg_c2v(ic,j) = c2v(ic,j);
      }
    }
    // Need an int**, even if only have one element type
    conn[0] = tg_c2v.getData();
  }
  else {
    // Need an int**, even if only have one element type
    conn[0] = c2v.getData();
  }

  tg->registerGridData(gridID,nVerts,xv.getData(),iblank.data(),nwall,nover,iwall.data(),
                       iover.data(),ntypes,nodesPerCell,&nEles,&conn[0]);

  // Create a new communicator to communicate across grids
  MPI_Comm_split(MPI_COMM_WORLD, gridRank, gridID, &interComm);
  MPI_Comm_split(MPI_COMM_WORLD, gridID, params->rank, &gridComm);
#endif
}

void geo::updateOversetConnectivity(void)
{
#ifndef _NO_MPI
  // Pre-process the grids
  tg->profile();

  // Have TIOGA perform the nodal overset connectivity (set nodal iblanks)
  tg->performConnectivity();

  // Only needed for debugging, really
  if (params->writeIBLANK)
    writeOversetConnectivity();

  // Now use new nodal iblanks to set cell and face iblanks
  setCellFaceIblanks();
#endif
}

void geo::setCellFaceIblanks()
{
  // Use the TIOGA-supplied nodal iblank values, set iblank values for all cells and faces

  // Only needed for moving grids: List of current hole cells
  holeCells.clear();
  blankCells.clear();
  unblankCells.clear();
  for (int ic=0; ic<nEles; ic++)
    if (iblankCell[ic] == HOLE)
      holeCells.insert(ic);

  iblankCell.assign(nEles,NORMAL);

  // First, blank all cells which contain a hole node

  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2nv[ic]; j++) {
      int iv = c2v(ic,j);
      if (iblank[iv] == HOLE) {
        iblankCell[ic] = HOLE;

        // Only needed for moving grids: Cells which must be removed from solver
        if (holeCells.count(ic))
          blankCells.insert(ic);

        break;
      }
    }
  }

  // Only needed for moving grids: Get cells which  must be 'un-blanked'
  for (auto &ic:holeCells)
    if (iblankCell[ic] == NORMAL)
      unblankCells.insert(ic);

  // Only needed for moving grids: List of current hole faces
  holeFaces.clear();
  blankFaces.clear();
  unblankFaces.clear();
  unblankOFaces.clear();
  for (int ic=0; ic<nFaces; ic++)
    if (iblankFace[ic] == HOLE)
      holeFaces.insert(ic);

  iblankFace.assign(nFaces,NORMAL);

  // Next, get the new overset faces & set all hole faces

  for (int ic=0; ic<nEles; ic++) {
    if (iblankCell[ic] == HOLE) {
      for (int j=0; j<c2nf[ic]; j++) {
        int ff = c2f(ic,j);
        if (iblankFace[ff] == NORMAL) {
          // If not set yet, assume fringe (overset), but set to hole if any hole nodes
          iblankFace[ff] = FRINGE;
          for (int k=0; k<f2nv[ff]; k++) {
            int iv = f2v(ff,k);
            if (iblank[iv] == HOLE) {
              iblankFace[ff] = HOLE;
              break;
            }
          }
        }
      }
    }
  }

  // Only needed for moving grids: Get faces which  must be 'un-blanked'
  for (auto &ic:holeFaces) {
    if (iblankFace[ic] == NORMAL) {
      unblankFaces.insert(ic);
    }
    else if (iblankFace[ic] == FRINGE) {
      unblankOFaces.insert(ic);
    }
  }

}

void geo::writeOversetConnectivity(void)
{
#ifndef _NO_MPI
  // Write out only the mesh with IBLANK info (no solution data)
  tg->writeData(0,NULL,0);
#endif
}

void geo::matchOversetDonors(vector<ele> &eles, vector<superMesh> &donors)
{
#ifndef _NO_MPI

#endif
}

vector<ele> geo::setupUnblankCells()
{
  vector<ele> unblankEles(unblankCells.size());

  int i=0;
  for (auto &ic:unblankCells) {
    unblankEles[i].ID = ic;
    // ...continue like with setupElesFaces()
    i++;
  }

  return unblankEles;
}
