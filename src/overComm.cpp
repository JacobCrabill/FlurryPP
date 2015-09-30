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

void overComm::setIblanks2D(matrix<double>& xv, matrix<int> &overFaces, matrix<int>& wallFaces, vector<int>& iblank)
{
  int nVerts = xv.getDim0();
  int nWallFaces = wallFaces.getDim0();
  int nOverFaces = overFaces.getDim0();

  Array<double,3> wallNodes(nWallFaces,2,2);
  for (int i=0; i<nWallFaces; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
        wallNodes(i,j,k) = xv(wallFaces(i,j),k);

  Array<double,3> overNodes(nOverFaces,2,2);
  for (int i=0; i<nOverFaces; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
        overNodes(i,j,k) = xv(overFaces(i,j),k);

  vector<int> nWallFace_rank, nOverFace_rank;
  vector<double> wallNodes_rank, overNodes_rank;

  gatherData(nWallFaces, 4, wallNodes.getData(), nWallFace_rank, wallNodes_rank);
  gatherData(nOverFaces, 4, overNodes.getData(), nOverFace_rank, overNodes_rank);

  /* --- Get bounding box of wall and overset nodes from all proccesses --- */

  point minPtW, maxPtW, minPtO, maxPtO;
  int offsetW = 0, offsetO = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) {
      offsetW += nWallFace_rank[p-1];
      offsetO += nOverFace_rank[p-1];
    }

    if (gridIdList[p] == gridID) continue;

    for (int i=0; i<nWallFace_rank[p]; i++) {
      for (int j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
          minPtW[k] = std::min(minPtW[k],wallNodes_rank[4*(offsetW+i)+2*j+k]);
          maxPtW[k] = std::max(maxPtW[k],wallNodes_rank[4*(offsetW+i)+2*j+k]);
        }
      }
    }

    for (int i=0; i<nOverFace_rank[p]; i++) {
      for (int j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
          minPtO[k] = std::min(minPtO[k],overNodes_rank[4*(offsetO+i)+2*j+k]);
          maxPtO[k] = std::max(maxPtO[k],overNodes_rank[4*(offsetO+i)+2*j+k]);
        }
      }
    }
  }

  /* --- Use winding-number method to find hole points given wall faces --- */

  iblank.assign(nVerts,NORMAL);

  double eps = 1e-3;
  double tol = 1e-10;
  for (int i=0; i<nVerts; i++) {
    point pt;
    pt.x = xv(i,0);
    pt.y = xv(i,1);

    int offsetW = 0;
    int offsetO = 0;
    double windW = 0;
    double windO = 0;
    for (int p=0; p<nproc; p++) {
      if (p>0) {
        offsetW += nWallFace_rank[p-1];
        offsetO += nOverFace_rank[p-1];
      }

      if (gridIdList[p] == gridID) continue;

      if ( (pt.x>minPtW.x-tol) && (pt.y>minPtW.y-tol) && (pt.x<maxPtW.x+tol) && (pt.y<maxPtW.y+tol) ) {
        // Point lies within bounding box of wall boundary, so calculate winding number
        for (int i=0; i<nWallFace_rank[p]; i++) {
          point pt1 = point(&wallNodes_rank[4*(offsetW+i)+0],2);
          point pt2 = point(&wallNodes_rank[4*(offsetW+i)+2],2);

          Vec3 dx1 = pt1 - pt;  dx1 /= dx1.norm();
          Vec3 dx2 = pt2 - pt;  dx2 /= dx2.norm();
          double dot = min(max(dx1*dx2,-1.),1.);
          windW += std::abs(std::acos(dot));
        }
      }

      if ( (pt.x>minPtO.x+tol) && (pt.y>minPtO.y+tol) && (pt.x<maxPtO.x-tol) && (pt.y<maxPtO.y-tol) ) {
        // Point lies within bounding box of overset boundary, minus tolerance, so calculate winding number
        for (int i=0; i<nOverFace_rank[p]; i++) {
          point pt1 = point(&overNodes_rank[4*(offsetO+i)+0],2);
          point pt2 = point(&overNodes_rank[4*(offsetO+i)+2],2);

          Vec3 dx1 = pt1 - pt;  dx1 /= dx1.norm();
          Vec3 dx2 = pt2 - pt;  dx2 /= dx2.norm();
          double dot = min(max(dx1*dx2,-1.),1.);
          windO += std::abs(std::acos(dot));
        }
      }
    }
    if (std::abs(windW)+eps >= 2.*PI)
      iblank[i] = HOLE;

    if (std::abs(windO)+eps >= 2.*PI && iblank[i]!=HOLE)
      iblank[i] = FRINGE;
  }
}

set<int> overComm::findCellDonors2D(vector<shared_ptr<ele>> &eles, const vector<double> &targetBox)
{
  // Find all eles which overlap with targetBox
  set<int> hitCells;
  for (auto &e:eles) {
    auto box = e->getBoundingBox();
    bool hit = true;
    for (int dim=0; dim<2; dim++) {
      hit = hit && (targetBox[dim+3] >= box[dim]);
      hit = hit && (targetBox[dim] <= box[dim+3]);
    }
    if (hit)
      hitCells.insert(e->ID);
  }

  return hitCells;
}

void overComm::matchOversetPoints3D(vector<shared_ptr<ele>> &eles, vector<shared_ptr<overFace>> &overFaces, const vector<int> &eleMap)
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

  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nPts_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    foundPts[p].resize(0);
    foundEles[p].resize(0);
    foundLocs[p].resize(0);
    for (int i=0; i<nPts_rank[p]; i++) {
      // Get requested interpolation point
      point pt = point(&interpPtsPhys[3*(offset+i)]);
      int ic = tg->findPointDonor(&interpPtsPhys[3*(offset+i)]);
      if (ic>=0 && eleMap[ic]>=0) {
        int ie = eleMap[ic];
        point refLoc;
        bool isInEle = eles[ie]->getRefLocNelderMeade(pt,refLoc);

        if (!isInEle) FatalError("Unable to match fringe point!");

        foundPts[p].push_back(i);
        foundEles[p].push_back(ic); // Local ele id for this grid
        foundLocs[p].push_back(refLoc);
      }
    }
  }

  /* ---- Prepare for Data Communication ---- */

  // Get the number of matched points for each grid
  nPtsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (gridIdList[p] == gridID) continue;
    nPtsSend[p] = foundPts[p].size();
  }
#endif
}

void overComm::matchOversetPoints2D(vector<shared_ptr<ele>> &eles, vector<shared_ptr<overFace>> &overFaces, const point &minPt, const point &maxPt)
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
  double tol = 1e-6;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nPts_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    foundPts[p].resize(0);
    foundEles[p].resize(0);
    foundLocs[p].resize(0);
    for (int i=0; i<nPts_rank[p]; i++) {
      // Get requested interpolation point
      point pt = point(&interpPtsPhys[3*(offset+i)]);

      // First, check that point even lies within bounding box of grid
      if ( (pt.x<minPt.x-tol) || (pt.y<minPt.y-tol) ||
           (pt.x>maxPt.x+tol) || (pt.y>maxPt.y+tol) )
        continue;

      // Check for containment in all eles on this rank of this grid [no ADT currently for 2D meshes]
      for (auto &e:eles) {
        point refLoc;
        bool isInEle = e->getRefLocNelderMeade(pt,refLoc);

        if (isInEle) {
          foundPts[p].push_back(i);
          foundEles[p].push_back(e->ID); // Local ele id for this grid
          foundLocs[p].push_back(refLoc);
          break;
        }
      }
    }
  }

  /* ---- Prepare for Data Communication ---- */

  // Get the number of matched points for each grid
  nPtsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (gridIdList[p] == gridID) continue;
    nPtsSend[p] = foundPts[p].size();
  }
#endif
}


void overComm::matchUnblankCells(vector<shared_ptr<ele>> &eles, map<int,map<int,oper>> &opers, set<int> &unblankCells, vector<int> &eleMap, int quadOrder)
{
#ifndef _NO_MPI

  /* ---- Send Unblanked-Cell Nodes to All Grids ---- */

  nUnblanks = unblankCells.size();

  int nDims = params->nDims;
  int nv = (nDims==2) ? 4 : 8;

  vector<int> ubCells;
  Array<double,3> ubCellNodes(nUnblanks,nv,nDims);
  int i = 0;
  for (auto &ic:unblankCells) {
    int ie = eleMap[ic];
    ubCells.push_back(ie);
    // Constraining this to just linear hexahedrons/quadrilaterals for the time being
    for (int j=0; j<nv; j++) {
      for (int k=0; k<nDims; k++) {
        ubCellNodes(i,j,k) = eles[ie]->nodesRK[j][k];
      }
    }
    i++;
  }

  /* ---- Gather all cell bounding-box data on each grid ---- */

  int stride = nv*nDims;
  vector<int> nCells_rank;
  vector<double> ubNodes_rank;
  gatherData(nUnblanks, stride, ubCellNodes.getData(), nCells_rank, ubNodes_rank);

  if (getSum(nCells_rank)==0) return;

  /* ---- Check Every Unblanked Cell for Donor Cells on This Grid ---- */

  foundCells.resize(nproc);
  foundCellDonors.resize(nproc);
  foundCellNDonors.resize(nproc);
  donors.resize(0);
  int offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nCells_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    foundCells[p].resize(0);
    foundCellDonors[p].setup(0,0);
    foundCellNDonors[p].resize(0);

    for (int i=0; i<nCells_rank[p]; i++) {
      // Get requested cell's bounding box [min & max extents]
      vector<double> targetBox = {1e15, 1e15, 1e15, -1e15, -1e15, -1e15};
      vector<point> targetNodes;
      for (int j=0; j<nv; j++) {
        point pt = point(&ubNodes_rank[(offset+i)*stride+nDims*j],nDims);
        targetNodes.push_back(pt);
        for (int dim=0; dim<3; dim++) {
          targetBox[dim]   = min(pt[dim],targetBox[dim]);
          targetBox[dim+3] = max(pt[dim],targetBox[dim+3]);
        }
      }

      // Find all possible donors using Tioga's ADT search (3D) or my brute-force search (2D)
      set<int> cellIDs;
      if (nDims == 2)
        cellIDs = findCellDonors2D(eles,targetBox);
      else
        cellIDs = tg->findCellDonors(targetBox.data());

      if (cellIDs.size() > 0) {
        vector<int> donorsIDs;
        foundCellNDonors[p].push_back(cellIDs.size());
        foundCells[p].push_back(i);
        for (auto &ic:cellIDs)
          donorsIDs.push_back(ic);

        // Setup the donor cells [on this grid] for the unblanked cell [on other grid]
        foundCellDonors[p].insertRowUnsized(donorsIDs);

        Array2D<point> donorPts;
        for (auto &ic:donorsIDs)
          donorPts.insertRow(eles[eleMap[ic]]->nodesRK);

        superMesh mesh(targetNodes,donorPts,quadOrder,nDims,rank,donors.size());
        donors.push_back(mesh);
      }
    }
  }

  /* --- Setup & Exchange Quadrature-Point Data --- */

  // Now that we have the local superMesh for each target, setup points for
  // use with Galerkin projection
  for (auto &mesh:donors)
    mesh.setupQuadrature();

  // Get the locations of the quadrature points for each target cell, and
  // the reference location of the points within the donor cells
  vector<matrix<double>> qpts(nproc), qptsD_ref(nproc), donorBasis(nproc), donorU(nproc);
  vector<vector<int>> targetID(nproc), donorID(nproc);
  int nSpts = (params->order+1)*(params->order+1);
  if (nDims == 3) nSpts *= (params->order+1);
  offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += foundCells[p-1].size();
    for (int i=0; i<foundCells[p].size(); i++) {
      vector<int> parents_tmp;
      matrix<double> qpts_tmp;
      donors[offset+i].getQpts(qpts_tmp,parents_tmp);

      for (int j=0; j<parents_tmp.size(); j++) {
        qpts[p].insertRow(qpts_tmp[j],INSERT_AT_END,3);
        donorID[p].push_back(foundCellDonors[p](i,parents_tmp[j]));
        targetID[p].push_back(foundCells[p][i]);
        int ic = eleMap[donorID[p].back()];
        point refLoc;
        bool isInEle = eles[ic]->getRefLocNelderMeade(point(qpts_tmp[j],nDims),refLoc);

        if (!isInEle) {
          cout << "qpt: " << qpts_tmp(j,0) << ", " << qpts_tmp(j,1) << endl;
          auto box = eles[ic]->getBoundingBox();
          cout << "ele box: " << box[0] << ", " << box[1] << "; " << box[3] << ", " << box[4] << endl;
          FatalError("Quadrature Point Reference Location not found in ele!");
        }
        qptsD_ref[p].insertRow({refLoc.x,refLoc.x,refLoc.z});

        vector<double> basisTmp;
        opers[eles[ic]->eType][eles[ic]->order].getBasisValues(refLoc,basisTmp);
        donorBasis[p].insertRow(basisTmp);
      }

      for (int id=0; id<foundCellNDonors[p][i]; id++) {
        int ic = eleMap[foundCellDonors[p](i,id)];
        for (int spt=0; spt<nSpts; spt++)
          donorU[p].insertRow(eles[ic]->U_spts[spt],INSERT_AT_END,nFields);
      }
    }
  }

  // Exchange superMesh quadrature points among the grids
  nQptsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (gridIdList[p] == gridID) continue;
    nQptsSend[p] = qpts[p].getDim0();
  }
  setupNPieces(nQptsSend,nQptsRecv);

  vector<vector<int>> unblankID(nproc);
  vector<matrix<double>> qpts_recv(nproc);
  sendRecvData(nQptsSend,nQptsRecv,targetID,unblankID,qpts,qpts_recv,3);

  /* --- Perform the Galerkin Projection --- */

  // Using target cell data and nodes, get basis-function and solution values
  // at all quadrature points
  vector<matrix<double>> sendBasis(nproc); // Basis function values to be sent back to donor grid
  offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nQptsRecv[p-1];
    if (gridIdList[p] == gridID) continue;

    for (int i=0; i<nQptsRecv[p]; i++) {
      int ii = unblankID[p][i];
      int ie = ubCells[ii];

      point refLoc;
      point pos = point(qpts_recv[p][i],nDims);
      bool isInEle = eles[ie]->getRefLocNelderMeade(pos,refLoc);
      if (!isInEle) FatalError("Quadrature Point Reference Location not found in ele!");

      qpts_recv[p](i,0) = refLoc.x;
      qpts_recv[p](i,1) = refLoc.y;
      qpts_recv[p](i,2) = refLoc.z;

      vector<double> basisTmp;
      opers[eles[ie]->eType][eles[ie]->order].getBasisValues(refLoc,basisTmp);

      sendBasis[p].insertRow(basisTmp);
    }
  }

  vector<matrix<double>> targetBasis(nproc);
  stride = nSpts;
  sendRecvData(nQptsRecv,nQptsSend,sendBasis,targetBasis,stride);

  // Perform integration for each target cell
  vector<matrix<double>> LHS(nproc), RHS(nproc);
  offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += foundCells[p-1].size();
    int offsetQ = 0;
    int offsetD = 0;
    for (int i=0; i<foundCells[p].size(); i++) {
      if (i>0) {
        offsetQ += donors[offset+i-1].getNQpts();
        offsetD += nSpts*foundCellNDonors[p][i-1];
      }

      matrix<double> rhs(nSpts,nFields);
      matrix<double> lhs;
      int nQpts = donors[offset+i].getNQpts();

      for (int ispt=0; ispt<nSpts; ispt++) {
        matrix<double> basisTgtDnr(nQpts,nSpts), basisTgtTgt(nQpts,nSpts);
        for (int qpt=0; qpt<nQpts; qpt++) {
          for (int jspt=0; jspt<nSpts; jspt++) {
            basisTgtDnr(qpt,jspt) = targetBasis[p](offsetQ+qpt,ispt) * donorBasis[p](offsetQ+qpt,jspt);
            basisTgtTgt(qpt,jspt) = targetBasis[p](offsetQ+qpt,ispt) * targetBasis[p](offsetQ+qpt,jspt);
          }
        }

        auto massMatTTRow = donors[offset+i].integrate(basisTgtTgt);
        lhs.insertRow(massMatTTRow);

        auto massMatTDRow = donors[offset+i].integrateByDonor(basisTgtDnr);

        for (int id=0; id<foundCellNDonors[p][i]; id++)
          for (int jspt=0; jspt<nSpts; jspt++)
            for (int k=0; k<nFields; k++)
              rhs(ispt,k) += massMatTDRow(id,jspt) * donorU[p](offsetD+id*nSpts+jspt,k);
      }
      LHS[p].appendRows(lhs);
      RHS[p].appendRows(rhs);
    }
  }

  // Send/recv the final target-cell data
  nCellsSend.resize(nproc);
  for (int p=0; p<nproc; p++) {
    nCellsSend[p] = foundCells[p].size();
  }
  setupNPieces(nCellsSend,nCellsRecv);

  int strideR = nSpts*nFields;
  int strideL = nSpts*nSpts;
  vector<matrix<double>> tmpUbLHS(nproc), tmpUbRHS(nproc);
  vector<vector<int>> recvInds;
  sendRecvData(nCellsSend,nCellsRecv,foundCells,recvInds,LHS,tmpUbLHS,strideL);
  sendRecvData(nCellsSend,nCellsRecv,RHS,tmpUbRHS,strideR);

  // Add contributions from superMeshes on each rank
  vector<matrix<double>> ubLHS(nUnblanks), ubRHS(nUnblanks);
  for (auto &lhs:ubLHS) lhs.setup(nSpts,nSpts);
  for (auto &rhs:ubRHS) rhs.setup(nSpts,nFields);

  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    for (int i=0; i<nCellsRecv[p]; i++) {
      for (int j=0; j<nSpts; j++) {
        for (int k=0; k<nSpts; k++) {
          ubLHS[recvInds[p][i]](j,k) += tmpUbLHS[p](i,j*nSpts+k);
        }
        for (int k=0; k<nFields; k++) {
          ubRHS[recvInds[p][i]](j,k) += tmpUbRHS[p](i,j*nFields+k);
        }
      }
    }
  }

  // Apply the new values to the unblank ele objects
  for (int i=0; i<nUnblanks; i++) {
    int ic = ubCells[i];
    eles[ic]->U_spts.initializeToZero();
    auto unblankU = solveCholesky(ubLHS[i],ubRHS[i]);
    for (int spt=0; spt<nSpts; spt++)
      for (int k=0; k<nFields; k++)
        eles[ic]->U_spts(spt,k) += unblankU(spt,k);
  }
#endif
}

void overComm::exchangeOversetData(vector<shared_ptr<ele>> &eles, map<int, map<int,oper> > &opers, vector<int> &eleMap)
{
#ifndef _NO_MPI

  U_out.resize(nproc);
  for (int p=0; p<nproc; p++) {
    U_out[p].setup(foundPts[p].size(),nFields);
    if (gridIdList[p] == gridID) continue;
    for (int i=0; i<foundPts[p].size(); i++) {
      point refPos = foundLocs[p][i];
      int ic = eleMap[foundEles[p][i]];
      if (ic<0 || ic>eles.size()) {
        cout << "!!!! ic = " << ic << " !!!!" << endl;
        cout << "rank " << params->rank << ", cell " << foundEles[p][i] << endl;
        FatalError("bad value of ic!");
      }
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
    cout << "rank " << params->rank << ", # Unmatched Points = " << nOverPts - getSum(nPtsRecv) << " out of " << nOverPts << endl;
    FatalError("Unmatched points remaining!");
  }

  recvPts.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    recvPts[p].resize(nPtsRecv[p]);
  }

  sendRecvData(nPtsSend,nPtsRecv,foundPts,recvPts,U_out,U_in,nFields,true);
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
#ifndef _NO_MPI
  nPiecesOut.assign(nproc,0);

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
#endif
}

template<typename T>
void overComm::sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<vector<int>> &recvInds,
                            vector<matrix<T>> &sendVals, matrix<T> &recvVals, int stride, bool matchInds=false)
{
  // Do basic send/receive, keeping receive values in destination arrays from each rank
  vector<matrix<T>> tmpRecvVals;
  if (matchInds)
    sendRecvData(nPiecesSend,nPiecesRecv,sendInds,recvInds,sendVals,tmpRecvVals,stride);
  else
    sendRecvData(nPiecesSend,nPiecesRecv,sendVals,tmpRecvVals,stride);

  // Rearrange data into final storage matrix
  recvVals.setup(getSum(nPiecesRecv),stride);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    for (int i=0; i<nPiecesRecv[p]; i++) {
      for (int k=0; k<stride; k++) {
        recvVals(recvInds[p][i],k) = tmpRecvVals[p](i,k);
      }
    }
  }
}

template<typename T>
void overComm::sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<matrix<T>> &sendVals, vector<matrix<T>> &recvVals, int stride)
{
#ifndef _NO_MPI
  MPI_Datatype T_TYPE = getMpiDatatype<T>();

  MPI_Status status;

  recvVals.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    recvVals[p].setup(nPiecesRecv[p],stride);
  }

  vector<MPI_Request> ValsRecvs(nproc);
  vector<MPI_Request> ValsSends(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    if (nPiecesRecv[p]>0)
      MPI_Irecv(recvVals[p].getData(),nPiecesRecv[p]*stride,T_TYPE,p,p,MPI_COMM_WORLD,&ValsRecvs[p]);
    if (nPiecesSend[p]>0)
      MPI_Isend(sendVals[p].getData(),nPiecesSend[p]*stride,T_TYPE,p,rank,MPI_COMM_WORLD,&ValsSends[p]);
  }

  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    if (nPiecesRecv[p]>0)
      MPI_Wait(&ValsRecvs[p], &status);
    if (nPiecesSend[p]>0)
      MPI_Wait(&ValsSends[p], &status);
  }
#endif
}

template<typename T>
void overComm::sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<vector<int>> &recvInds,
                            vector<matrix<T>> &sendVals, vector<matrix<T>> &recvVals, int stride)
{
#ifndef _NO_MPI
  MPI_Datatype T_TYPE = getMpiDatatype<T>();

  MPI_Status status;

  recvInds.resize(nproc);
  recvVals.resize(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    recvInds[p].resize(nPiecesRecv[p]);
    recvVals[p].setup(nPiecesRecv[p],stride);
  }

  vector<MPI_Request> IndsRecvs(nproc);
  vector<MPI_Request> IndsSends(nproc);
  vector<MPI_Request> ValsRecvs(nproc);
  vector<MPI_Request> ValsSends(nproc);
  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    if (nPiecesRecv[p]>0) {
      MPI_Irecv(recvInds[p].data(),nPiecesRecv[p],MPI_INT,p,p,MPI_COMM_WORLD,&IndsRecvs[p]);
      MPI_Irecv(recvVals[p].getData(),nPiecesRecv[p]*stride,T_TYPE,p,p,MPI_COMM_WORLD,&ValsRecvs[p]);
    }
    if (nPiecesSend[p]>0) {
      MPI_Isend(sendInds[p].data(),nPiecesSend[p],MPI_INT,p,rank,MPI_COMM_WORLD,&IndsSends[p]);
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
#endif
}
