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

void overComm::matchOversetPoints(vector<shared_ptr<ele>> &eles, vector<shared_ptr<overFace>> &overFaces, matrix<int> &c2ac,
                                  const vector<int> &eleMap, const point &centroid, const point &extents)
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

  /* ---- Check Every Fringe Point for Donor Cell on This Grid ---- */

  foundPts.resize(nproc);
  foundEles.resize(nproc);
  foundLocs.resize(nproc);
  int offset = 0;

  /*double eps = 1e-10;
  point minPt, maxPt;
  minPt.x = centroid.x - extents.x/2. - eps;
  minPt.y = centroid.y - extents.y/2. - eps;
  minPt.z = centroid.z - extents.z/2. - eps;

  maxPt.x = centroid.x + extents.x/2. + eps;
  maxPt.y = centroid.y + extents.y/2. + eps;
  maxPt.z = centroid.z + extents.z/2. + eps;*/

  /*for (int p=0; p<nproc; p++) {
    if (p>0) offset += nPts_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    int currCell = 0;  // Global ele ID
    int ic = eleMap[currCell];        // Corresponding index in 'eles' vector
    for (int i=0; i<nPts_rank[p]; i++) {
      // Get current interpolation point
      point pt = point(&interpPtsPhys[3*(offset+i)]);

      // First, check that point even lies within bounding box of grid
      if ( (pt.x<minPt.x) || (pt.y<minPt.y) || (pt.z<minPt.z) ||
           (pt.x>maxPt.x) || (pt.y>maxPt.y) || (pt.z>maxPt.z) )
        continue;

      // Find the next valid, non-blanked element to start with
      while (ic<0) {
        currCell++;
        if (currCell>=c2ac.getDim0()) currCell = 0;
        ic = eleMap[currCell];
      }

      point refLoc;
      bool isInEle = eles[ic]->getRefLocNelderMeade(pt,refLoc);

      set<int> triedCells;
      triedCells.insert(ic);

      // Setup first layer of cells to search
      set<int> nextCells;
      for (int j=0; j<c2ac.getDim1(); j++) {
        int ic0 = c2ac(currCell,j);
        if (ic0>0 && eleMap[ic0]>0)
          nextCells.insert(ic0);
      }

      // Use c2c to march towards interpolation point
      while (!isInEle && nextCells.size()>0) {
        double minDist = 1e15;
        int nextCell = -1;
        for (auto &icNext:nextCells) {
          int icn = eleMap[icNext];
          isInEle = eles[icn]->getRefLocNelderMeade(pt,refLoc);
          triedCells.insert(icn);

          if (isInEle) {
            currCell = icNext;
            ic = icn;
            break;
          }
          else {
            auto box = eles[icn]->getBoundingBox();
            point xc_next = point(&box[0]);
            double dist = getDist(pt,xc_next);
            if (dist < minDist) {
              minDist = dist;
              nextCell = icNext;
              currCell = icNext;
            }
          }
        }

        if (isInEle || nextCell<0) break;

        // Setup the next layer of cells to search based on the best cell from the previous layer
        nextCells.clear();
        for (int j=0; j<c2ac.getDim1(); j++) {
          int icn = c2ac(nextCell,j);
          if (icn>0 && eleMap[icn]>0 && !triedCells.count(eleMap[icn]))
            nextCells.insert(icn);
        }
      }

      if (isInEle) {
        foundPts[p].push_back(i);
        foundEles[p].push_back(ic); // Local ele id for this grid
        foundLocs[p].push_back(refLoc);
      }
    }
  }*/

  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nPts_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    for (int i=0; i<nPts_rank[p]; i++) {
      // Get requested interpolation point
      point pt = point(&interpPtsPhys[3*(offset+i)]);
      int ic = tg->findPointDonor(&interpPtsPhys[3*(offset+i)]);
      if (ic>=0 && eleMap[ic]>=0) {
        ic = eleMap[ic];
        point refLoc;
        eles[ic]->getRefLocNelderMeade(pt,refLoc);
        foundPts[p].push_back(i);
        foundEles[p].push_back(ic); // Local ele id for this grid
        foundLocs[p].push_back(refLoc);
      }

      /*// First, check that point even lies within bounding box of grid
      if ( (pt.x<minPt.x) || (pt.y<minPt.y) || (pt.z<minPt.z) ||
           (pt.x>maxPt.x) || (pt.y>maxPt.y) || (pt.z>maxPt.z) )
        continue;

      // Check for containment in all eles on this rank of this grid
      int ic = 0;
      for (auto &e:eles) {
        // !! THIS IS HORRIBLY INEFFICIENT !! - should use c2c to march from an
        // initial guess to the neighboring cell nearest the point in question
        point refLoc;
        bool isInEle = e->getRefLocNelderMeade(pt,refLoc);

        if (isInEle) {
          foundPts[p].push_back(i);
          foundEles[p].push_back(ic); // Local ele id for this grid
          foundLocs[p].push_back(refLoc);
          break;
        }
        ic++;
      }*/
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
MPI_Barrier(MPI_COMM_WORLD);
  U_in.setup(nOverPts,nFields);
#endif
}

void overComm::matchUnblankCells(vector<shared_ptr<ele>> &eles, set<int> &unblankCells, matrix<int> &c2c, vector<int> &eleMap, int quadOrder)
{
#ifndef _NO_MPI

  /* ---- Send Unblanked-Cell Nodes to All Grids ---- */

  nUnblanks = unblankCells.size();

  Array<double,3> ubCellNodes(nUnblanks,8,3);
  int i = 0;
  for (auto &ic:unblankCells) {
    int ie = eleMap[ic];
    // Constraining this to just linear hexahedrons for the time being
    for (int j=0; j<8; j++) {
      for (int k=0; k<3; k++) {
        ubCellNodes(i,j,k) = eles[ie]->nodesRK[j][k];
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
      // Get requested cell's bounding box [min & max extents]
      vector<double> targetBox = {1e15, 1e15, 1e15, -1e15, -1e15, -1e15};
      vector<point> targetNodes;
      for (int j=0; j<8; j++) {
        point pt = point(&ubNodes_rank[(offset+i)*stride+3*j]);
        targetNodes.push_back(pt);
        for (int dim=0; dim<3; dim++) {
          targetBox[dim]   = min(pt[dim],targetBox[dim]);
          targetBox[dim+3] = max(pt[dim],targetBox[dim+3]);
        }
      }

      // Find all possible donors using Tioga's ADT search
      set<int> cellIDs = tg->findCellDonors(targetBox.data());

      if (cellIDs.size() > 0) {
        vector<int> donorsIDs;
        int ind = foundCells.size();
        foundCellNDonors[p].push_back(cellIDs.size());
        foundCells[p].push_back(i);
        for (auto &ic:cellIDs)
          donorsIDs.push_back(ic);

        // Setup the donor cells [on this grid] for the unblanked cell [on other grid]
        foundCellDonors[p].insertRowUnsized(donorsIDs);

        Array2D<point> donorPts;
        for (auto &ic:donorsIDs) {
          donorPts.insertRow(eles[eleMap[ic]]->nodesRK);
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

void overComm::exchangeOversetData(vector<shared_ptr<ele>> &eles, map<int, map<int,oper> > &opers)
{
#ifndef _NO_MPI

  U_out.resize(nproc);
  for (int p=0; p<nproc; p++) {
    U_out[p].setup(foundPts[p].size(),nFields);
    if (gridIdList[p] == gridID) continue;
    for (int i=0; i<foundPts[p].size(); i++) {
      point refPos = foundLocs[p][i];
      int ic = foundEles[p][i];
      opers[eles[ic]->eType][eles[ic]->order].interpolateToPoint(eles[ic]->U_spts, U_out[p][i], refPos);
    }
  }

  /* ---- Send/Receive the the interpolated data across grids using interComm ---- */
  nPtsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (gridIdList[p] == gridID) continue;
    nPtsSend[p] = foundPts[p].size();
  }

  setupNPieces(nPtsSend,nPtsRecv);

  if (nOverPts > getSum(nPtsRecv)) {
    cout << "rank " << rank << ", # Unmatched Points = " << nOverPts - getSum(nPtsRecv) << " out of " << nOverPts << endl;
    FatalError("Unmatched points remaining!");
  }

  recvPts.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    recvPts[p].resize(nPtsRecv[p]);
  }

  sendRecvData(nPtsSend,nPtsRecv,foundPts,recvPts,U_out,U_in,nFields,1);

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
