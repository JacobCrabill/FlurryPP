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
  tg = new tioga();

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

  iblankCell.assign(nEles,NORMAL);

  // First, blank all cells which contain a hole node
  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2nv[ic]; j++) {
      int iv = c2v(ic,j);
      if (iblank[iv] == HOLE) {
        iblankCell[ic] = HOLE;
        break;
      }
    }
  }

  // Next, get the new overset faces & set all hole faces
  iblankFace.assign(nFaces,NORMAL);
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

}

void geo::writeOversetConnectivity(void)
{
#ifndef _NO_MPI
  // Write out only the mesh with IBLANK info (no solution data)
  tg->writeData(0,NULL,0);
#endif
}

void geo::matchOversetPoints(vector<ele> &eles, vector<shared_ptr<overFace>> &overFacesVec)
{
#ifndef _NO_MPI
  // Get all of the fringe points on this grid
  overPts.resize(0);
  for (auto &oface: overFacesVec) {
    oface->fptOffset = overPts.size();
    auto pts = oface->getPosFpts();
    overPts.insert(overPts.end(),pts.begin(),pts.end());
  }

  overPtsPhys = createMatrix(overPts);
  int nInterpPts = overPts.size();

  /* ---- Gather all interpolation point data on each grid ---- */

  // Create a new communicator to communicate across grids
  MPI_Comm_split(MPI_COMM_WORLD, gridRank, gridID, &interComm);
  MPI_Comm_split(MPI_COMM_WORLD, gridID, params->rank, &gridComm);

  // Get the number of interp points for each process on this grid
  nPts_rank.resize(nprocPerGrid);
  MPI_Allgather(&nInterpPts, 1, MPI_INT, nPts_rank.data(), 1, MPI_INT, gridComm);

  // Accumulate all interpolation points on this grid into single matrix
  int nPtsGrid = 0;
  vector<int> recvCnts(nprocPerGrid);
  vector<int> recvDisp(nprocPerGrid);
  for (int i=0; i<nprocPerGrid; i++) {
    nPtsGrid += nPts_rank[i];
    recvCnts[i] = nPts_rank[i]*3;
    if (i>0)
      recvDisp[i] = recvDisp[i-1] + recvCnts[i-1];
  }

  if (gridRank == 0)
    cout << "Geo: Grid " << gridID << ": # of Overset Points = " << nPtsGrid << endl;

  matrix<double> overPts_rank(nPtsGrid,3);
  MPI_Allgatherv(overPtsPhys.getData(), overPtsPhys.getSize(), MPI_DOUBLE, overPts_rank.getData(), recvCnts.data(), recvDisp.data(), MPI_DOUBLE, gridComm);

  // Send this rank's data to all other grids
  vector<int> nPts_grid(nGrids);
  MPI_Allgather(&nInterpPts, 1, MPI_INT, nPts_grid.data(), 1, MPI_INT, interComm);

  // Reduce across each grid
  vector<int> nPts_tmp = nPts_grid;
  MPI_Allreduce(nPts_tmp.data(), nPts_grid.data(), nGrids, MPI_INT, MPI_SUM, gridComm);

  // Setup storage for all interpolation points in simulation
  int nPtsTotal = getSum(nPts_grid);
  interpPtsPhys.setup(nPtsTotal,3);

  recvCnts.resize(nGrids);
  recvDisp.assign(nGrids,0);
  for (int i=0; i<nGrids; i++) {
    recvCnts[i] = nPts_grid[i]*3;
    if (i>0)
      recvDisp[i] = recvDisp[i-1] + recvCnts[i-1];
  }
  // Send this grid's points to each grid along gridRank
  MPI_Allgatherv(overPts_rank.getData(), overPts_rank.getSize(), MPI_DOUBLE, interpPtsPhys.getData(), recvCnts.data(), recvDisp.data(), MPI_DOUBLE, interComm);

  /* ---- Check Every Fringe Point for Donor Cell on This Grid ---- */

  foundPts.resize(nGrids);
  foundEles.resize(nGrids);
  foundLocs.resize(nGrids);
  int offset = 0;
  for (int g=0; g<nGrids; g++) {
    if (g>0) offset += nPts_grid[g-1];

    if (g == gridID) continue;

    for (int i=0; i<nPts_grid[g]; i++) {
      // Get requested interpolation point
      point pt = point(interpPtsPhys[offset+i]);
      // Check for containment in all eles on this rank of this grid
      for (auto &e:eles) {
        /* !! THIS IS HORRIBLY INEFFICIENT !! - should use c2c to march from an
         * initial guess to the neighboring cell nearest the point in question */
        point refLoc;
        bool isInEle = e.getRefLocNelderMeade(pt,refLoc);

        if (isInEle) {
          foundPts[g].push_back(i);
          foundEles[g].push_back(e.ID); // Local ele id for this grid
          foundLocs[g].push_back(refLoc);
          break;
        }
      }
    }
  }

  /* ---- Send/Receive the the donor info across grids ---- */

  // Send/Receive the number of matched points on each grid
  nPtsRecv.resize(nGrids);
  nPtsSend.resize(nGrids);
  vector<MPI_Request> reqs(nGrids);
  vector<MPI_Request> sends(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    nPtsSend[g] = foundPts[g].size();
    MPI_Irecv(&nPtsRecv[g], 1, MPI_INT, g, gridID, interComm, &reqs[g]);
    MPI_Isend(&nPtsSend[g], 1, MPI_INT, g, g, interComm, &sends[g]);
  }

  MPI_Status status;
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    MPI_Wait(&reqs[g],&status);
    MPI_Wait(&sends[g],&status);
  }
#endif
}


void geo::exchangeOversetData(vector<matrix<double>> &U_ipts, matrix<double> &U_opts)
{
#ifndef _NO_MPI
  /* ---- Send/Receive the the interpolated data across grids using interComm ---- */

  // Send/Receive the number of matched points on each grid
  vector<MPI_Request> recvsPts(nGrids);
  vector<MPI_Request> sendsPts(nGrids);
  vector<MPI_Request> recvsVals(nGrids);
  vector<MPI_Request> sendsVals(nGrids);

  // Send/Receive the accompanying point IDs which were matched
  vector<vector<int>> recvPts(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    recvPts[g].resize(nPtsRecv[g]);
  }

  vector<matrix<double>> recvU(nGrids);
  for (uint g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    recvU[g].setup(nPtsRecv[g],nFields);
  }

  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    if (nPtsRecv[g] > 0) {
      MPI_Irecv(recvPts[g].data(),  nPtsRecv[g],        MPI_INT,    g, gridID, interComm, &recvsPts[g]);
      MPI_Irecv(recvU[g].getData(), recvU[g].getSize(), MPI_DOUBLE, g, gridID, interComm, &recvsVals[g]);
    }
    if (nPtsSend[g] > 0) {
      MPI_Isend(foundPts[g].data(),  nPtsSend[g],         MPI_INT,    g, g, interComm, &sendsPts[g]);
      MPI_Isend(U_ipts[g].getData(), U_ipts[g].getSize(), MPI_DOUBLE, g, g, interComm, &sendsVals[g]);
    }
  }

  MPI_Status status;
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    if (nPtsRecv[g] > 0) {
      MPI_Wait(&recvsPts[g],&status);
      MPI_Wait(&recvsVals[g],&status);
    }
    if (nPtsSend[g] > 0) {
      MPI_Wait(&sendsPts[g],&status);
      MPI_Wait(&sendsVals[g],&status);
    }
  }

  /* ---- Distribute the the interpolated data within each grid ---- */

  // Now that each grid has the matched point donor info distributed among its
  // processes, send to proper gridRank
  vector<int> nGPtsSend(nprocPerGrid);
  vector<int> nGPtsRecv(nprocPerGrid);
  vector<vector<int>> gridPtIdsSend(nprocPerGrid);
  vector<vector<int>> gridPtIdsRecv(nprocPerGrid);
  vector<matrix<double>> gridURecv(nprocPerGrid);
  vector<matrix<double>> gridUSend(nprocPerGrid);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    for (int i=0; i<nPtsRecv[g]; i++) {
      int iv = recvPts[g][i];
      int nv = nPts_rank[0];
      int p = 0;
      int offset = 0;
      while (nv<=iv && p<nprocPerGrid) {
        offset += nPts_rank[p];
        p++;
        nv += nPts_rank[p];
      }
      iv -= offset; // Make iv be the process-local point ID
      gridPtIdsSend[p].push_back(iv);
      auto U = recvU[g].getRow(i);
      gridUSend[p].insertRow(U);
      nGPtsSend[p]++;
    }
  }

  vector<MPI_Request> reqsGrid(nprocPerGrid);
  vector<MPI_Request> sendsGrid(nprocPerGrid);
  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Irecv(&nGPtsRecv[p], 1, MPI_INT, p, gridRank, gridComm, &reqsGrid[p]);
    MPI_Isend(&nGPtsSend[p], 1, MPI_INT, p, p,        gridComm, &sendsGrid[p]);
  }

  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Wait(&reqsGrid[p],&status);
    MPI_Wait(&sendsGrid[p],&status);
  }

  // Send/Receive the accompanying point IDs and gridIDs which were matched
  for (int p=0; p<nprocPerGrid; p++) {
    gridPtIdsRecv[p].resize(nGPtsRecv[p]);
    gridURecv[p].setup(nGPtsRecv[p],nFields);
  }

  vector<MPI_Request> reqsPt(nprocPerGrid);
  vector<MPI_Request> sendsPt(nprocPerGrid);
  vector<MPI_Request> reqsU(nprocPerGrid);
  vector<MPI_Request> sendsU(nprocPerGrid);
  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Irecv(gridPtIdsRecv[p].data(), nGPtsRecv[p],           MPI_INT,    p, gridRank, gridComm, &reqsPt[p]);
    MPI_Irecv(gridURecv[p].getData(),  gridURecv[p].getSize(), MPI_DOUBLE, p, gridRank, gridComm, &reqsU[p]);
    MPI_Isend(gridPtIdsSend[p].data(), nGPtsSend[p],           MPI_INT,    p, p, gridComm, &sendsPt[p]);
    MPI_Isend(gridUSend[p].getData(),  gridUSend[p].getSize(), MPI_DOUBLE, p, p, gridComm, &sendsU[p]);
  }

  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Wait(&reqsPt[p],&status);
    MPI_Wait(&reqsU[p],&status);
    MPI_Wait(&sendsPt[p],&status);
    MPI_Wait(&sendsU[p],&status);
  }

  /* ---- Gather all data together into output data storage ---- */

  for (int p=0; p<nprocPerGrid; p++) {
    for (int i=0; i<nGPtsRecv[p]; i++) {
      int iv = gridPtIdsRecv[p][i];
      for (int k=0; k<nFields; k++) {
        U_opts(iv,k) = gridURecv[p](i,k);
      }
    }
  }
#endif
}
