/*!
 * \file overComm.cpp
 * \brief Source file for overComm (Overset Communicator) class
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

void overComm::matchOversetPoints(vector<ele> &eles, vector<shared_ptr<overFace>> &overFaces)
{
#ifndef _NO_MPI
  /* ---- Gather all interpolation point data on each grid ---- */

  // Get all of the fringe points on this grid
  overPts.setup(0,0);
  for (auto &oface: overFaces) {
    oface->OComm = this;
    oface->fptOffset = overPts.getDim0();
    auto pts = oface->getPosFpts();
    for (auto &pt:pts)
      overPts.insertRow({pt.x, pt.y, pt.z});
  }

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

  /* ---- Prepare for Data Communication ---- */

  // Get the number of matched points for each grid
  nPtsSend.resize(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    nPtsSend[g] = foundPts[g].size();
  }

  U_in.setup(nOverPts,nFields);
#endif
}

void overComm::matchOversetUnblanks(vector<ele> &eles, vector<shared_ptr<overFace>> &overFaces, set<int> &unblankCells, set<int>& blankedOFaces,
                                    set<int> &unblankOFaces, vector<int> &eleMap, vector<int> &faceMap, int quadOrder)
{
#ifndef _NO_MPI

  /* --- Find amount of data to be sent around --- */

  nUnblanks = unblankCells.size();

  /* ---- Send Unblanked-Cell Nodes to All Grids ---- */

  Array<double,3> ubCellNodes(nUnblanks,8,3);
  int i = 0;
  for (auto &ic:unblankCells) {
    int ie = eleMap[ic];
    // Constraining this to just linear hexahedrons for the time being
    for (int j=0; j<8; j++) {
      for (int k=0; k<3; k++) {
        ubCellNodes(i,j,k) = eles[ie].nodesRK[0][j][k];
      }
    }
    i++;
  }

  /* ---- Gather all cell bounding-box data on each grid ---- */

  // Gather bounding boxes for all unblanked cells on this rank
  matrix<double> bBoxes;
  for (auto &ic:unblankCells) {
    int ie = eleMap[ic];
    bBoxes.insertRow(eles[ie].getBoundingBox());
  }

  int stride = 24;
  vector<int> nCells_rank, nCells_grid;
  vector<double> ubNodes_grid;
  gatherData(nUnblanks, stride, ubCellNodes.getData(), nCells_rank, nCells_grid, ubNodes_grid);

  /* ---- Check Every Unblanked Cell for Donor Cells on This Grid ---- */

  foundCells.resize(nGrids);
  foundCellDonors.resize(nGrids);
  foundCellNDonors.resize(nGrids);
  int offset = 0;
  for (int g=0; g<nGrids; g++) {
    if (g>0) offset += nCells_grid[g-1];

    if (g == gridID) continue;

    for (int i=0; i<nCells_grid[g]; i++) {
      // Get requested cell's bounding box
      vector<point> targetNodes(8);
      for (int j=0; j<8; j++) {
        targetNodes[j] = point(&ubNodes_grid[(offset+i)*stride+3*j]);
      }
      point cent, dx;
      getBoundingBox(targetNodes,cent,dx);
      //point cent = point(&bBoxes_grid[offset+i]);
      //point dx = point(&bBoxes_grid[offset+i+3]);

      // Check for overlap in all eles on this rank of this grid
      bool found = false;
      int ind = -1;
      vector<int> donorsIDs;
      Array2D<point> donorPts;
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
          donorsIDs.push_back(e.ID); // Local ele id for this grid
        }
      }
      if (found) {
        // Setup the donor cells [on this grid] for the unblanked cell [on other grid]
        foundCellDonors[g].insertRowUnsized(donorsIDs);
        for (int k=0; k<donorsIDs.size(); k++) {
          donorPts.insertRow(eles[donorsIDs[k]].nodesRK[0]);
        }
        superMesh mesh(targetNodes,donorPts,quadOrder);
        donors.push_back(mesh);
      }
    }
  }

  /* ---- Now that we have the local superMesh for each target, setup points for
   *      use with Galerkin projection ---- */

  for (auto &mesh:donors)
    mesh.setupQuadrature();

  // Next Up: Match the new overset-face points in a separate set of data arrays,
  // remove the 'blanked' overset face points from the original data arrays, then
  // interlace the new data into the original data
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
  nPtsSend.resize(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    nPtsSend[g] = foundPts[g].size();
  }

  distributeData(nPtsSend,nPtsRecv,foundPts,nPts_rank,U_out,U_in,nFields);
#endif
}

#ifndef _NO_MPI
template<typename T>
MPI_Datatype overComm::getMpiDatatype(void)
{
  FatalError("MPI Datatype not implemented here");
}

template<> MPI_Datatype overComm::getMpiDatatype<int>(void){ return MPI_INT; }

template<> MPI_Datatype overComm::getMpiDatatype<double>(void){ return MPI_DOUBLE; }

template<> MPI_Datatype overComm::getMpiDatatype<float>(void){ return MPI_FLOAT; }

template<> MPI_Datatype overComm::getMpiDatatype<unsigned int>(void){ return MPI_UNSIGNED; }

template<> MPI_Datatype overComm::getMpiDatatype<long>(void){ return MPI_LONG; }

template<> MPI_Datatype overComm::getMpiDatatype<short>(void){ return MPI_SHORT; }

template<> MPI_Datatype overComm::getMpiDatatype<char>(void){ return MPI_CHAR; }
#endif

template<typename T>
void overComm::gatherData(int nPieces, int stride, T *values, vector<int>& nPieces_rank, vector<int> &nPieces_grid, vector<T> &values_all)
{
#ifndef _NO_MPI

  /* ---- Gather all overset data pieces on each grid ---- */

  MPI_Datatype T_TYPE = getMpiDatatype<T>();

  // Get the number of data pieces for each process on this grid
  nPieces_rank.resize(nprocPerGrid);
  MPI_Allgather(&nPieces, 1, MPI_INT, nPieces_rank.data(), 1, MPI_INT, gridComm);

  // Accumulate all data pieces on this grid into single matrix
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
  MPI_Allgatherv(values, nPieces*stride, T_TYPE, values_rank.data(), recvCnts.data(), recvDisp.data(), T_TYPE, gridComm);

  // Send this rank's data to all other grids
  nPieces_grid.resize(nGrids);
  MPI_Allgather(&nPieces, 1, MPI_INT, nPieces_grid.data(), 1, MPI_INT, interComm);

  // Reduce across each grid
  vector<int> nPieces_tmp = nPieces_grid;
  MPI_Allreduce(nPieces_tmp.data(), nPieces_grid.data(), nGrids, MPI_INT, MPI_SUM, gridComm);

  // Setup storage for all data pieces across all grids/ranks
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
  // Send this grid's data to each grid along gridRank
  MPI_Allgatherv(values_rank.data(), values_rank.size(), T_TYPE, values_all.data(), recvCnts.data(), recvDisp.data(), T_TYPE, interComm);
#endif
}

template<typename T>
void overComm::distributeData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<int> &nPieces_rank, vector<matrix<T>> &values_send, matrix<T> &values_recv, int stride)
{
#ifndef _NO_MPI
  MPI_Datatype T_TYPE = getMpiDatatype<T>();

  MPI_Status status;

  // Send/Receive the number of data pieces on each grid
  nPiecesRecv.resize(nGrids);
  vector<MPI_Request> reqs(nGrids);
  vector<MPI_Request> sends(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    MPI_Irecv(&nPiecesRecv[g], 1, MPI_INT, g, gridID, interComm, &reqs[g]);
    MPI_Isend(&nPiecesSend[g], 1, MPI_INT, g, g, interComm, &sends[g]);
  }

  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    MPI_Wait(&reqs[g],&status);
    MPI_Wait(&sends[g],&status);
  }

  // Send/Receive the accompanying piece destination indices

  vector<MPI_Request> recvsInds(nGrids);
  vector<MPI_Request> sendsInds(nGrids);
  vector<MPI_Request> recvsVals(nGrids);
  vector<MPI_Request> sendsVals(nGrids);

  vector<vector<int>> recvInds(nGrids);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    recvInds[g].resize(nPiecesRecv[g]);
  }

  vector<matrix<T>> recvU(nGrids);
  for (uint g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    recvU[g].setup(nPiecesRecv[g],stride);
  }

  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    if (nPiecesRecv[g] > 0) {
      MPI_Irecv(recvInds[g].data(), nPiecesRecv[g],     MPI_INT, g, gridID, interComm, &recvsInds[g]);
      MPI_Irecv(recvU[g].getData(), recvU[g].getSize(), T_TYPE,  g, gridID, interComm, &recvsVals[g]);
    }
    if (nPiecesSend[g] > 0) {
      MPI_Isend(sendInds[g].data(),       nPiecesSend[g],           MPI_INT, g, g, interComm, &sendsInds[g]);
      MPI_Isend(values_send[g].getData(), values_send[g].getSize(), T_TYPE,  g, g, interComm, &sendsVals[g]);
    }
  }

  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    if (nPiecesRecv[g] > 0) {
      MPI_Wait(&recvsInds[g],&status);
      MPI_Wait(&recvsVals[g],&status);
    }
    if (nPiecesSend[g] > 0) {
      MPI_Wait(&sendsInds[g],&status);
      MPI_Wait(&sendsVals[g],&status);
    }
  }

  /* ---- Distribute the data within each grid ---- */

  // Now that each grid has the data and destination info distributed among its
  // processes, send to proper gridRank
  vector<int> nGPcsSend(nprocPerGrid);
  vector<int> nGPcsRecv(nprocPerGrid);
  vector<vector<int>> gridIndsSend(nprocPerGrid);
  vector<vector<int>> gridIndsRecv(nprocPerGrid);
  vector<matrix<T>> gridValsRecv(nprocPerGrid);
  vector<matrix<T>> gridValsSend(nprocPerGrid);
  for (int g=0; g<nGrids; g++) {
    if (g == gridID) continue;
    for (int i=0; i<nPiecesRecv[g]; i++) {
      int ind = recvInds[g][i];
      int nv = nPieces_rank[0];
      int p = 0;
      int offset = 0;
      while (nv<=ind && p<nprocPerGrid) {
        offset += nPieces_rank[p];
        p++;
        nv += nPieces_rank[p];
      }
      ind -= offset; // Make ind be the process-local value ID
      gridIndsSend[p].push_back(ind);
      auto U = recvU[g].getRow(i);
      gridValsSend[p].insertRow(U);
      nGPcsSend[p]++;
    }
  }

  vector<MPI_Request> recvsGrid(nprocPerGrid);
  vector<MPI_Request> sendsGrid(nprocPerGrid);
  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Irecv(&nGPcsRecv[p], 1, MPI_INT, p, gridRank, gridComm, &recvsGrid[p]);
    MPI_Isend(&nGPcsSend[p], 1, MPI_INT, p, p,        gridComm, &sendsGrid[p]);
  }

  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Wait(&recvsGrid[p],&status);
    MPI_Wait(&sendsGrid[p],&status);
  }

  // Send/Receive the accompanying point IDs and gridIDs which were matched
  for (int p=0; p<nprocPerGrid; p++) {
    gridIndsRecv[p].resize(nGPcsRecv[p]);
    gridValsRecv[p].setup(nGPcsRecv[p],stride);
  }

  recvsInds.resize(nprocPerGrid);
  sendsInds.resize(nprocPerGrid);
  recvsVals.resize(nprocPerGrid);
  sendsVals.resize(nprocPerGrid);
  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Irecv(gridIndsRecv[p].data(),    nGPcsRecv[p],              MPI_INT, p, gridRank, gridComm, &recvsInds[p]);
    MPI_Irecv(gridValsRecv[p].getData(), gridValsRecv[p].getSize(), T_TYPE,  p, gridRank, gridComm, &recvsVals[p]);
    MPI_Isend(gridIndsSend[p].data(),    nGPcsSend[p],              MPI_INT, p, p, gridComm, &sendsInds[p]);
    MPI_Isend(gridValsSend[p].getData(), gridValsSend[p].getSize(), T_TYPE,  p, p, gridComm, &sendsVals[p]);
  }

  for (int p=0; p<nprocPerGrid; p++) {
    MPI_Wait(&recvsInds[p],&status);
    MPI_Wait(&recvsVals[p],&status);
    MPI_Wait(&sendsInds[p],&status);
    MPI_Wait(&sendsVals[p],&status);
  }

  /* ---- Gather all data together into output data storage
   *      (Note that values_recv should be pre-sized correctly) ---- */

  for (int p=0; p<nprocPerGrid; p++) {
    for (int i=0; i<nGPcsRecv[p]; i++) {
      int iv = gridIndsRecv[p][i];
      for (int k=0; k<stride; k++) {
        values_recv(iv,k) = gridValsRecv[p](i,k);
      }
    }
  }
#endif
}

