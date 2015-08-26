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

void overComm::setup(input* _params, int _nGrids, int _gridID, int _gridRank, int _nprocPerGrid, vector<int>& _gridIdList)
{
  params = _params;

  nGrids = _nGrids;
  gridID = _gridID;
  gridRank = _gridRank;
  nprocPerGrid = _nprocPerGrid;
  gridIdList = _gridIdList;

  rank = params->rank;
  nproc = params->nproc;
  nFields = params->nFields;

#ifndef _NO_MPI
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

  vector<double> interpPtsPhys;

  gatherData(nOverPts, 3, overPts.getData(), nPts_rank, interpPtsPhys);

  //cout << "Grid,Rank = " << gridID << "," << rank << ": nOverPts = " << nOverPts << endl;

  /* ---- Check Every Fringe Point for Donor Cell on This Grid ---- */

  foundPts.resize(nproc);
  foundEles.resize(nproc);
  foundLocs.resize(nproc);
  int offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nPts_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    for (int i=0; i<nPts_rank[p]; i++) {
      // Get requested interpolation point
      point pt = point(&interpPtsPhys[3*(offset+i)]);
      //point pt = point(interpPtsPhys[offset+i]);
      // Check for containment in all eles on this rank of this grid
      int ic = 0;
      for (auto &e:eles) {
        /* !! THIS IS HORRIBLY INEFFICIENT !! - should use c2c to march from an
         * initial guess to the neighboring cell nearest the point in question */
        point refLoc;
        bool isInEle = e.getRefLocNelderMeade(pt,refLoc);

        if (isInEle) {
          foundPts[p].push_back(i);
          foundEles[p].push_back(ic); // Local ele id for this grid
          foundLocs[p].push_back(refLoc);
          break;
        }
        ic++;
      }
    }
    //cout << "Grid,Rank = " << gridID << "," << rank << ": nFoundPts for rank" << p << " = " << foundPts[p].size() << endl;
  }

  /* ---- Prepare for Data Communication ---- */

  // Get the number of matched points for each grid
  nPtsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (gridIdList[p] == gridID) continue;
    nPtsSend[p] = foundPts[p].size();
  }

  U_in.setup(nOverPts,nFields);
#endif
}

void overComm::matchUnblankCells(vector<ele> &eles, set<int> &unblankCells, vector<int> &eleMap, int quadOrder)
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

  int stride = 24;
  vector<int> nCells_rank;
  vector<double> ubNodes_rank;
  gatherData(nUnblanks, stride, ubCellNodes.getData(), nCells_rank, ubNodes_rank);

  /* ---- Check Every Unblanked Cell for Donor Cells on This Grid ---- */

  foundCells.resize(nproc);
  foundCellDonors.resize(nproc);
  foundCellNDonors.resize(nproc);
  int offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nCells_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    for (int i=0; i<nCells_rank[p]; i++) {
      // Get requested cell's bounding box
      vector<point> targetNodes(8);
      for (int j=0; j<8; j++) {
        targetNodes[j] = point(&ubNodes_rank[(offset+i)*stride+3*j]);
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
            foundCellNDonors[p].push_back(0);
            foundCells[p].push_back(i);
            found = true;
          }
          foundCellNDonors[p][ind]++;
          donorsIDs.push_back(e.ID); // Local ele id for this grid
        }
      }
      if (found) {
        // Setup the donor cells [on this grid] for the unblanked cell [on other grid]
        foundCellDonors[p].insertRowUnsized(donorsIDs);
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

  U_out.resize(nproc);
  for (int p=0; p<nproc; p++) {
    U_out[p].setup(foundPts[p].size(),nFields);
    if (gridIdList[p] == gridID) continue;
    for (int i=0; i<foundPts[p].size(); i++) {
      point refPos = foundLocs[p][i];
      int ic = foundEles[p][i];
      opers[eles[ic].eType][eles[ic].order].interpolateToPoint(eles[ic].U_spts, U_out[p][i], refPos);
    }
  }

  /* ---- Send/Receive the the interpolated data across grids using interComm ---- */
  nPtsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (gridIdList[p] == gridID) continue;
    nPtsSend[p] = foundPts[p].size();
  }

  setupNPieces(nPtsSend,nPtsRecv);

//  int sum=0;
//  for (int p=0; p<nproc; p++) {
//    if (gridIdList[p] == gridID) continue;
//    sum += nPtsRecv[p];
//  }
//  if (sum != nOverPts) {
//  //if (getSum(nPtsRecv) != nOverPts) {
//    cout << "rank = " << rank << ": diff = " << nOverPts-getSum(nPtsRecv) << ": ";
//    FatalError("Did not get all overset points matched!");
//  }

  recvPts.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    recvPts[p].resize(nPtsRecv[p]);
  }

  sendRecvData(nPtsSend,nPtsRecv,foundPts,recvPts,U_out,U_in,nFields,1);

//  if (U_in.checkNan()) cout << "rank " << rank << ": NaN in U_in" << endl;

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
void overComm::gatherData(int nPieces, int stride, T *values, vector<int>& nPieces_rank, vector<T> &values_all)
{
#ifndef _NO_MPI

  /* ---- Gather all overset data pieces on each grid ---- */

  MPI_Datatype T_TYPE = getMpiDatatype<T>();

  // Get the number of data pieces for each process on this grid
  nPieces_rank.resize(nproc);
  MPI_Allgather(&nPieces, 1, MPI_INT, nPieces_rank.data(), 1, MPI_INT, MPI_COMM_WORLD);

  int nPiecesTotal = 0;
  vector<int> recvCnts(nproc);
  vector<int> recvDisp(nproc);
  for (int p=0; p<nproc; p++) {
    nPiecesTotal += nPieces_rank[p];
    recvCnts[p] = nPieces_rank[p]*stride;
    if (p>0)
      recvDisp[p] = recvDisp[p-1] + recvCnts[p-1];
  }

  values_all.resize(nPiecesTotal*stride);
  MPI_Allgatherv(values, nPieces*stride, T_TYPE, values_all.data(), recvCnts.data(), recvDisp.data(), T_TYPE, MPI_COMM_WORLD);
#endif
}

void overComm::setupNPieces(vector<int> &nPiecesIn, vector<int> &nPiecesOut)
{
  nPiecesOut.resize(nproc);

  vector<MPI_Request> sends(nproc);
  vector<MPI_Request> recvs(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    MPI_Irecv(&nPiecesOut[p], 1, MPI_INT, p, p,    MPI_COMM_WORLD, &recvs[p]);
    MPI_Isend(&nPiecesIn[p],  1, MPI_INT, p, rank, MPI_COMM_WORLD, &sends[p]);
  }

  MPI_Status status;
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    MPI_Wait(&recvs[p], &status);
    MPI_Wait(&sends[p], &status);
  }
}

template<typename T>
void overComm::sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<vector<int>> &recvInds,
                            vector<matrix<T>> &sendVals, matrix<T> &recvVals, int stride, int matchInds=0)
{
#ifndef _NO_MPI
  MPI_Datatype T_TYPE = getMpiDatatype<T>();

  MPI_Status status;

  vector<matrix<T>> tmpRecvVals(nproc);

  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    tmpRecvVals[p].setup(nPtsRecv[p],stride);
  }

  vector<MPI_Request> IndsRecvs(nproc);
  vector<MPI_Request> IndsSends(nproc);
  vector<MPI_Request> ValsRecvs(nproc);
  vector<MPI_Request> ValsSends(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    if (nPiecesRecv[p]>0) {
      if (matchInds) {
        // If needed, receive the local array indices which the data is matched with
        MPI_Irecv(recvInds[p].data(),nPiecesRecv[p],MPI_INT,p,p,MPI_COMM_WORLD,&IndsRecvs[p]);
      }
      MPI_Irecv(tmpRecvVals[p].getData(),nPiecesRecv[p]*stride,T_TYPE,p,p,MPI_COMM_WORLD,&ValsRecvs[p]);
    }
    if (nPiecesSend[p]>0) {
      if (matchInds) {
        // If needed, send the destination indices which the data is matched with
        MPI_Isend(sendInds[p].data(),nPiecesSend[p],MPI_INT,p,rank,MPI_COMM_WORLD,&IndsSends[p]);
      }
      MPI_Isend(sendVals[p].getData(),nPiecesSend[p]*stride,T_TYPE,p,rank,MPI_COMM_WORLD,&ValsSends[p]);
    }
  }

  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    if (nPiecesRecv[p]>0) {
      MPI_Wait(&IndsRecvs[p], &status);
      MPI_Wait(&ValsRecvs[p], &status);
    }
    if (nPiecesSend[p]>0) {
      MPI_Wait(&IndsSends[p], &status);
      MPI_Wait(&ValsSends[p], &status);
    }
  }

  // Rearrange data into final storage matrix  [NOTE: U_in must be pre-sized]
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    for (int i=0; i<nPtsRecv[p]; i++) {
      for (int k=0; k<stride; k++) {
        recvVals(recvInds[p][i],k) = tmpRecvVals[p](i,k);
      }
    }
  }
#endif
}
