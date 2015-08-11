/*!
 * \file overComm.hpp
 * \brief Header file for overComm class
 *
 * Handles communication of data across multiple MPI-partitioned overset grids
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

#include "overComm.hpp"

#include "global.hpp"

overComm::overComm()
{

}

void overComm::setup(input* _params, int _nGrids, int _gridID, int _gridRank, int _nprocPerGrid)
{
  params = _params;

  nGrids = _nGrids;
  gridID = _gridID;
  gridRank = _gridRank;
  nprocPerGrid = _nprocPerGrid;

  nFields = params->nFields;

#ifndef _NO_MPI
  MPI_Comm_split(MPI_COMM_WORLD, gridRank, gridID, &interComm);
  MPI_Comm_split(MPI_COMM_WORLD, gridID, params->rank, &gridComm);
#endif
}

void overComm::matchOversetPoints(vector<ele> &eles)
{
#ifndef _NO_MPI
  /* ---- Gather all interpolation point data on each grid ---- */

  nOverPts = overPts.getDim0();

  vector<int> nPts_grid(nGrids);
  vector<double> interpPtsPhys;

  gatherData(nOverPts, 3, overPts.getData(), nPts_rank, nPts_grid, interpPtsPhys);

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
      point pt = point(&interpPtsPhys[3*(offset+i)]);
      //point pt = point(interpPtsPhys[offset+i]);
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

  U_in.setup(nOverPts,nFields);
#endif
}

void overComm::matchOversetUnblanks(vector<ele> &eles, set<int> &unblankCells)
{
#ifndef _NO_MPI

  /* --- Find amount of data to be sent around --- */

  nUnblanks = unblankCells.size();

  /* ---- Gather all cell bounding-box data on each grid ---- */

  // Gather bounding boxes for all unblanked cells on this rank
  matrix<double> bBoxes;
  for (auto &e:eles) {
    if (unblankCells.count(e.ID)) {
      bBoxes.insertRow(e.getBoundingBox());
    }
  }

  vector<int> nCells_rank, nCells_grid;
  vector<double> bBoxes_grid;
  gatherData(nUnblanks, 1, bBoxes.getData(), nCells_rank, nCells_grid, bBoxes_grid);

  /* ---- Check Every Fringe Point for Donor Cell on This Grid ---- */

  foundCells.resize(nGrids);
  foundCellDonors.resize(nGrids);
  foundCellNDonors.resize(nGrids);
  int offset = 0;
  for (int g=0; g<nGrids; g++) {
    if (g>0) offset += nCells_grid[g-1];

    if (g == gridID) continue;

    for (int i=0; i<nCells_grid[g]; i++) {
      // Get requested cell's bounding box
      point cent = point(&bBoxes_grid[offset+i]);
      point dx = point(&bBoxes_grid[offset+i+3]);

      // Check for overlap in all eles on this rank of this grid
      bool found = false;
      int ind = -1;
      vector<int> donors;
      for (auto &e:eles) {
        /* !! THIS IS HORRIBLY INEFFICIENT !! - should use c2c to march from an
         * initial guess to the neighboring cell nearest the point in question */
/*     !! Sketch for using c2c: !!
        set<int> triedCells;
        int currCell = 0;
        for (int i=0; i<nCells_grid[g]; i++) {
          bool matched = false;
          while (!matched) {
            bool intersect = hitTest(box,currCell);
            triedCells.insert(currCell);
            if (intersect) {
              foundCellDonors[g].push_back(currCell);
              vector<int> possibleDonors = c2c.getRow(ic);
              // Explore around current cell to find other possible donors. If no intersection, we're too far away.
              while (possibleDonors.size()>0) {
                int nextC = possibleDonors.back();
                triedCells.insert(nextC);
                possibleDonors.pop_back();
                if (hitTest(bbox,eles[nextC])) {
                  foundCellDonors[g].push_back(nextC);
                  for (int j=0; j<c2nc[nextC]; j++) {
                    if (!triedCells.count(c2c(nextC,j))) {
                      possibleDonors.push_back(c2c(nextC,j));
                    }
                  }
                }
              }
            }
            else {
              // center of current bounding box
              point xc = point(boundingBoxex[offset+i]);
              point ele_xc1 = eles[currCell].getCentroid();
              double maxDot = 0;
              int nextCell = 0;
              for (int j=0; j<c2nc[ic]; j++) {
                point ele_xc2 = eles[c2c(ic,j)].getCentroid();
                Vec3 dx1 = xc - ele_xc1;
                Vec3 dx2 = ele_xc2 - ele_xc1;
                double dot = dx1*dx2;
                if (dot > maxDot) {
                  maxDot = dot;
                  nextCell = j;
                }
              }
              currCell = c2c(currCell,nextCell);
            }
          }
        }
*/

        auto eBox = e.getBoundingBox();

        bool intersect = (abs(cent.x-eBox[0])*2. < (dx.x + eBox[3])) &&
                         (abs(cent.y-eBox[1])*2. < (dx.y + eBox[4])) &&
                         (abs(cent.z-eBox[2])*2. < (dx.z + eBox[5]));

        if (intersect) {
          if (!found) {
            ind = foundCells.size();
            foundCellNDonors[g].push_back(0);
            foundCells[g].push_back(i);
            found = true;
          }
          foundCellNDonors[g][ind]++;
          donors.push_back(e.ID); // Local ele id for this grid
        }
      }
      foundCellDonors[g].insertRowUnsized(donors);
    }
  }

  /* ---- Send/Receive the the donor info across grids ---- */

  // Send/Receive the number of matched points on each grid
  nCellsRecv.resize(nGrids);
  nCellsSend.resize(nGrids);
  vector<MPI_Request> reqs(nGrids);
  vector<MPI_Request> sends(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    nCellsSend[g] = foundCells[g].size();
    MPI_Irecv(&nCellsRecv[g], 1, MPI_INT, g, gridID, interComm, &reqs[g]);
    MPI_Isend(&nCellsSend[g], 1, MPI_INT, g, g, interComm, &sends[g]);
  }

  MPI_Status status;
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    MPI_Wait(&reqs[g],&status);
    MPI_Wait(&sends[g],&status);
  }
#endif
}

template<typename T>
MPI_Datatype overComm::getMpiDatatype(void)
{
  FatalError("MPI Datatype not implemented here");
}

template<>
MPI_Datatype overComm::getMpiDatatype<int>(void)
{
  return MPI_INT;
}

template<>
MPI_Datatype overComm::getMpiDatatype<double>(void)
{
  return MPI_DOUBLE;
}

template<typename T>
void overComm::gatherData(int nPieces, int stride, T *values, vector<int>& nPieces_rank, vector<int> &nPieces_grid, vector<T> &values_all)
{
#ifndef _NO_MPI
  /* ---- Gather all interpolation point data on each grid ---- */

  MPI_Datatype T_Type = getMpiDatatype<T>();

  // Get the number of interp points for each process on this grid
  nPieces_rank.resize(nprocPerGrid);
  MPI_Allgather(&nPieces, 1, MPI_INT, nPieces_rank.data(), 1, MPI_INT, gridComm);

  // Accumulate all interpolation points on this grid into single matrix
  int nPiecesGrid = 0;
  vector<int> recvCnts(nprocPerGrid);
  vector<int> recvDisp(nprocPerGrid);
  for (int p=0; p<nprocPerGrid; p++) {
    nPiecesGrid += nPieces_rank[p];
    recvCnts[p] = nPieces_rank[p]*stride;
    if (p>0)
      recvDisp[p] = recvDisp[p-1] + recvCnts[p-1];
  }

  vector<T> values_rank(nPiecesGrid*stride);
  MPI_Allgatherv(values, nPieces*stride, T_Type, values_rank.data(), recvCnts.data(), recvDisp.data(), T_Type, gridComm);

  // Send this rank's data to all other grids
  nPieces_grid.resize(nGrids);
  MPI_Allgather(&nPieces, 1, MPI_INT, nPieces_grid.data(), 1, MPI_INT, interComm);

  // Reduce across each grid
  vector<int> nPieces_tmp = nPieces_grid;
  MPI_Allreduce(nPieces_tmp.data(), nPieces_grid.data(), nGrids, MPI_INT, MPI_SUM, gridComm);

  // Setup storage for all interpolation points in simulation
  int nPiecesTotal = getSum(nPieces_grid);

  values_all.resize(nPiecesTotal*stride);

  // If no processes have any data to send, return
  if (nPiecesTotal == 0) return;

  recvCnts.resize(nGrids);
  recvDisp.assign(nGrids,0);
  for (int g=0; g<nGrids; g++) {
    recvCnts[g] = nPieces_grid[g]*stride;
    if (g>0)
      recvDisp[g] = recvDisp[g-1] + recvCnts[g-1];
  }
  // Send this grid's points to each grid along gridRank
  MPI_Allgatherv(values_rank.data(), values_rank.size(), T_Type, values_all.data(), recvCnts.data(), recvDisp.data(), T_Type, interComm);
#endif
}

void overComm::exchangeOversetData(vector<ele> &eles, map<int, map<int,oper> > &opers)
{
#ifndef _NO_MPI

  U_out.resize(nGrids);
  for (int g=0; g<nGrids; g++) {
    U_out[g].setup(foundPts[g].size(),nFields);
    if (g == gridID) continue;
    for (int i=0; i<foundPts[g].size(); i++) {
      point refPos = foundLocs[g][i];
      int ic = foundEles[g][i];
      opers[eles[ic].eType][eles[ic].order].interpolateToPoint(eles[ic].U_spts, U_out[g][i], refPos);
    }
  }

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
      MPI_Isend(foundPts[g].data(), nPtsSend[g],        MPI_INT,    g, g, interComm, &sendsPts[g]);
      MPI_Isend(U_out[g].getData(), U_out[g].getSize(), MPI_DOUBLE, g, g, interComm, &sendsVals[g]);
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
        U_in(iv,k) = gridURecv[p](i,k);
      }
    }
  }
#endif
}

