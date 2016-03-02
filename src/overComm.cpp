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
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
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

#include "overComm.hpp"

#include "flux.hpp"
#include "global.hpp"

#ifndef _NO_MPI
template<typename T>
MPI_Datatype getMpiDatatype(void)
{
  FatalError("MPI Datatype not implemented here");
}

template<> MPI_Datatype getMpiDatatype<int>(void){ return MPI_INT; }

template<> MPI_Datatype getMpiDatatype<double>(void){ return MPI_DOUBLE; }

template<> MPI_Datatype getMpiDatatype<float>(void){ return MPI_FLOAT; }

template<> MPI_Datatype getMpiDatatype<unsigned int>(void){ return MPI_UNSIGNED; }

template<> MPI_Datatype getMpiDatatype<long>(void){ return MPI_LONG; }

template<> MPI_Datatype getMpiDatatype<short>(void){ return MPI_SHORT; }

template<> MPI_Datatype getMpiDatatype<char>(void){ return MPI_CHAR; }
#endif

overComm::overComm()
{

}

void overComm::setup(input* _params, int _nGrids, int _gridID, int _gridRank, int _nprocPerGrid, vector<int>& _gridIdList)
{
  params = _params;

  nDims = params->nDims;
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
    if (std::abs(windW)+eps >= 2.*pi)
      iblank[i] = HOLE;

    if (std::abs(windO)+eps >= 2*pi && iblank[i]!=HOLE)
      iblank[i] = FRINGE;
  }
}

unordered_set<int> overComm::findCellDonors2D(vector<shared_ptr<ele>> &eles, const vector<double> &targetBox)
{
  // Find all eles which overlap with targetBox
  unordered_set<int> hitCells;
  for (auto &e:eles) {
    auto box = e->getBoundingBox();
    bool hit = true;
    for (int dim=0; dim<2; dim++) {
      hit = hit && (targetBox[dim+2] >= box[dim]);
      hit = hit && (targetBox[dim] <= box[dim+3]);
    }
    if (hit)
      hitCells.insert(e->ID);
  }

  return hitCells;
}

void overComm::setupOverFacePoints(vector<shared_ptr<overFace>> &overFaces)
{
  // Get all of the fringe points on this grid
  overPts.setup(0,0);
  overNorm.setup(0,0);
  for (auto &oface: overFaces) {
    oface->OComm = this;
    oface->fptOffset = overPts.getDim0();
    auto pts = oface->getPosFpts();
    for (auto &pt:pts)
      overPts.insertRow({pt.x, pt.y, pt.z});
    if (params->oversetMethod==1) {
      auto norms = oface->getNormFpts();
      for (auto &vec:norms)
        overNorm.insertRow({vec.x,vec.y,vec.z});
    }
  }
}

void overComm::setupFringeCellPoints(vector<shared_ptr<ele>> &eles, const unordered_set<int> &fringeCells, const vector<int> &eleMap)
{
  overPts.setup(0,0);
  for (auto &ie:fringeCells) {
    int ic = eleMap[ie];
    eles[ic]->sptOffset = overPts.getDim0();
    auto pts = eles[ic]->getPosSpts();
    for (auto &pt:pts)
      overPts.insertRow({pt.x, pt.y, pt.z});
  }
}

void overComm::transferEleData(vector<shared_ptr<ele>> &eles, const unordered_set<int> &fringeCells, const vector<int> &eleMap)
{
  int nFields = params->nFields;

  for (auto &ie:fringeCells) {
    int ic = eleMap[ie];
    int Row = eles[ic]->sptOffset;
    for (int spt=0; spt<eles[ic]->nSpts; spt++) {
      for (int k=0; k<nFields; k++) {
        eles[ic]->U_spts(spt,k) = U_in(Row+spt,k);
      }
    }
  }
}

void overComm::matchOversetPoints(vector<shared_ptr<ele>> &eles, const vector<int> &eleMap, const point &minPt, const point &maxPt)
{
#ifndef _NO_MPI
  /* ---- Gather all interpolation point data on each grid ---- */

  nOverPts = overPts.getDim0();

  vector<double> interpPtsPhys;

  gatherData(nOverPts, 3, overPts.getData(), nPts_rank, interpPtsPhys);

  vector<double> interpNorms;
  if (params->oversetMethod==1)
    gatherData(nOverPts, 3, overNorm.getData(), nPts_rank, interpNorms);

  /* ---- Check Every Fringe Point for Donor Cell on This Grid ---- */

  // For use with ADT
  if (params->nDims==2) {
    if (eleList.size() != eleMap.size()) {
      eleList.resize(eleMap.size());
      for (int i=0; i<eleMap.size(); i++) eleList[i] = i;
    }
  }

  foundPts.resize(nproc);
  foundEles.resize(nproc);
  foundLocs.resize(nproc);
  if (params->oversetMethod==1)
    foundNorm.resize(nproc);
  int offset = 0;
  double tol = 1e-6;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nPts_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    foundPts[p].resize(0);
    foundEles[p].resize(0);
    foundLocs[p].resize(0);
    if (params->oversetMethod==1)
      foundNorm[p].resize(0);
    for (int i=0; i<nPts_rank[p]; i++) {
      // Get requested interpolation point
      point pt = point(&interpPtsPhys[3*(offset+i)]);

      if (params->nDims == 2) {
        // First, check that point even lies within bounding box of grid
        if ( (pt.x<minPt.x-tol) || (pt.y<minPt.y-tol) ||
             (pt.x>maxPt.x+tol) || (pt.y>maxPt.y+tol) )
          continue;

        // Use ADT to find all cells whose bounding box contains the point
        unordered_set<int> cellIDs;
        vector<double> targetBox = {pt.x,pt.y,pt.x,pt.y};
        adt->searchADT_box(eleList.data(),cellIDs,targetBox.data());
        for (auto &ic:cellIDs) {
          if (eleMap[ic]<0) continue;
          point refLoc;
          bool isInEle = eles[eleMap[ic]]->getRefLocNewton(pt,refLoc);

          if (isInEle) {
            foundPts[p].push_back(i);
            foundEles[p].push_back(ic); // Local ele id for this grid
            foundLocs[p].push_back(refLoc);
            if (params->oversetMethod==1)
              foundNorm[p].push_back(point(&interpNorms[3*(offset+i)]));
            break;
          }
        }
      }
      else {
        int ic = tg->findPointDonor(&interpPtsPhys[3*(offset+i)]);
        if (ic>=0 && eleMap[ic]>=0) {
          int ie = eleMap[ic];
          point refLoc;
          bool isInEle = eles[ie]->getRefLocNelderMead(pt,refLoc);

          if (!isInEle) FatalError("Unable to match fringe point!");

          foundPts[p].push_back(i);
          foundEles[p].push_back(ic); // Local ele id for this grid
          foundLocs[p].push_back(refLoc);
          if (params->oversetMethod==1)
            foundNorm[p].push_back(point(&interpNorms[3*(offset+i)]));
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


void overComm::matchUnblankCells(vector<shared_ptr<ele>> &eles, unordered_set<int>& unblankCells, vector<int> &eleMap, int quadOrder)
{
#ifndef _NO_MPI
  /* ---- Send Unblanked-Cell Nodes to All Grids ---- */

  nUnblanks = unblankCells.size();

  int nv = (nDims==2) ? 4 : 8;
  quadOrder = max(quadOrder,1);

  ubCells.resize(0);
  Array<double,3> ubCellNodes(nUnblanks,nv,nDims);
  int i = 0;
  for (auto &ic:unblankCells) {
    int ie = eleMap[ic];
    if (ie<0) {
      _print(params->rank,ic);
      FatalError("Unblank cell not in ele map...");
    }
    ubCells.push_back(ie);
    // Constraining this to just linear hexahedrons/quadrilaterals for the time being
    if (params->motion) {
      for (int j=0; j<nv; j++)
        for (int k=0; k<nDims; k++)
          ubCellNodes(i,j,k) = eles[ie]->nodesRK(j,k);
    } else {
      for (int j=0; j<nv; j++)
        for (int k=0; k<nDims; k++)
          ubCellNodes(i,j,k) = eles[ie]->nodes(j,k);
    }

    i++;
  }

  /* ---- Gather all cell bounding-box data on each grid ---- */

  int stride = nv*nDims;
  vector<int> nCells_rank;
  vector<double> ubNodes_rank;
  gatherData(nUnblanks, stride, ubCellNodes.getData(), nCells_rank, ubNodes_rank);

  nUnblanksTotal = getSum(nCells_rank);
  if (nUnblanksTotal==0) return;

  /* ---- Check Every Unblanked Cell for Donor Cells on This Grid ---- */

  if (eleList.size() != eleMap.size()) {
    eleList.resize(eleMap.size());
    for (int i=0; i<eleMap.size(); i++) eleList[i] = i;
  }

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
      vector<double> targetBox(nDims*2);
      for (int dim=0; dim<nDims; dim++) {
        targetBox[dim]       = 1e15;
        targetBox[dim+nDims] =-1e15;
      }
      vector<point> targetNodes;
      for (int j=0; j<nv; j++) {
        point pt = point(&ubNodes_rank[(offset+i)*stride+nDims*j],nDims);
        targetNodes.push_back(pt);
        for (int dim=0; dim<nDims; dim++) {
          targetBox[dim]       = min(pt[dim],targetBox[dim]);
          targetBox[dim+nDims] = max(pt[dim],targetBox[dim+nDims]);
        }
      }

      // Find all possible donors using Tioga's ADT search (3D) or my brute-force search (2D)
      unordered_set<int> cellIDs;
      if (nDims == 2) {
        adt->searchADT_box(eleList.data(),cellIDs,targetBox.data());
      }
      else {
        cellIDs = tg->findCellDonors(targetBox.data());
      }

      if (cellIDs.size() > 0) {
        vector<int> donorsIDs;
        foundCellNDonors[p].push_back(cellIDs.size());
        foundCells[p].push_back(i);
        for (auto &ic:cellIDs)
          donorsIDs.push_back(ic);

        // Setup the donor cells [on this grid] for the unblanked cell [on other grid]
        foundCellDonors[p].insertRowUnsized(donorsIDs);

        Array2D<point> donorPts;
        if (params->motion) {
          for (auto &ic:donorsIDs) {
            vector<point> tmpPts;
            for (uint npt = 0; npt < eles[eleMap[ic]]->nNodes; npt++)
              tmpPts.push_back(point(&eles[eleMap[ic]]->nodesRK(npt,0),nDims));
            donorPts.insertRow(tmpPts);
          }
        } else {
          for (auto &ic:donorsIDs) {
            vector<point> tmpPts;
            for (uint npt = 0; npt < eles[eleMap[ic]]->nNodes; npt++)
              tmpPts.push_back(point(&eles[eleMap[ic]]->nodes(npt,0),nDims));
            donorPts.insertRow(tmpPts);
          }
        }

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
#endif
}

void overComm::performGalerkinProjection(vector<shared_ptr<ele>> &eles, map<int, oper > &opers, vector<int> &eleMap, int order)
{
#ifndef _NO_MPI
  if (nUnblanksTotal == 0) return;

  // Get the locations of the quadrature points for each target cell, and
  // the reference location of the points within the donor cells
  qpts.resize(nproc);
  qptsD_ref.resize(nproc);
  donorBasis.resize(nproc);
  targetID.resize(nproc);
  donorID.resize(nproc);

  vector<matrix<double>> donorU(nproc);
  int nSpts = (order+1)*(order+1);
  if (nDims == 3) nSpts *= (order+1);
  int offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += foundCells[p-1].size();

    qpts[p].setup(0,0);
    qptsD_ref[p].setup(0,0);
    donorID[p].resize(0);
    targetID[p].resize(0);
    donorBasis[p].setup(0,0);

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
        bool isInEle;
        if (nDims==2) {
          isInEle = eles[ic]->getRefLocNewton(point(qpts_tmp[j],nDims),refLoc);
          if (!isInEle) {
            // Second try with slower (but more robust) algorithm
            isInEle = eles[ic]->getRefLocNelderMead(point(qpts_tmp[j],nDims),refLoc);
          }
        } else {
          isInEle = eles[ic]->getRefLocNelderMead(point(qpts_tmp[j],nDims),refLoc);
        }

        if (!isInEle)
          FatalError("Quadrature Point Reference Location not found in ele!");

        qptsD_ref[p].insertRow({refLoc.x,refLoc.x,refLoc.z});

        vector<double> basisTmp;
        opers[eles[ic]->order].getBasisValues(refLoc,basisTmp);
        donorBasis[p].insertRow(basisTmp);
      }

      for (int id=0; id<foundCellNDonors[p][i]; id++) {
        int ic = eleMap[foundCellDonors[p](i,id)];
        for (int spt=0; spt<nSpts; spt++)
          donorU[p].insertRow(&eles[ic]->U_spts(spt,0),INSERT_AT_END,nFields);
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
      bool isInEle = eles[ie]->getRefLocNelderMead(pos,refLoc);

      if (!isInEle) {
        cout << setprecision(16) << "qpt: " << pos << endl;
        cout << "Ele nodes:" << endl;
//        for (int n=0; n<4; n++)
//          _(eles[ie]->nodes[n]);
        FatalError("Quadrature Point Reference Location not found in ele!");
      }

      qpts_recv[p](i,0) = refLoc.x;
      qpts_recv[p](i,1) = refLoc.y;
      qpts_recv[p](i,2) = refLoc.z;

      vector<double> basisTmp;
      opers[eles[ie]->order].getBasisValues(refLoc,basisTmp);

      sendBasis[p].insertRow(basisTmp);
    }
  }

  vector<matrix<double>> targetBasis(nproc);
  int stride = nSpts;
  sendRecvData(nQptsRecv,nQptsSend,sendBasis,targetBasis,stride);

  // Perform integration for each target cell
  vector<matrix<double>> LHS(nproc), RHS(nproc);
  massMatTDRow.resize(donors.size()*nSpts);
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

        massMatTDRow[(offset+i)*nSpts+ispt] = donors[offset+i].integrateByDonor(basisTgtDnr);

        for (int id=0; id<foundCellNDonors[p][i]; id++)
          for (int jspt=0; jspt<nSpts; jspt++)
            for (int k=0; k<nFields; k++)
              rhs(ispt,k) += massMatTDRow[(offset+i)*nSpts+ispt](id,jspt) * donorU[p](offsetD+id*nSpts+jspt,k);
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
  sendRecvData(nCellsSend,nCellsRecv,foundCells,recvInds,LHS,tmpUbLHS,strideL);
  sendRecvData(nCellsSend,nCellsRecv,RHS,tmpUbRHS,strideR);

  // Add contributions from superMeshes on each rank
  vector<matrix<double>> ubRHS(nUnblanks);
  ubLHS.resize(nUnblanks);
  for (auto &lhs:ubLHS) {
    lhs.setup(nSpts,nSpts);
    lhs.initializeToZero();
  }
  for (auto &rhs:ubRHS) {
    rhs.setup(nSpts,nFields);
    rhs.initializeToZero();
  }

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
    for (int spt=0; spt<nSpts; spt++)
      for (int k=0; k<nFields; k++)
        eles[ic]->U_spts(spt,k) = 0;

    auto unblankU = solveCholesky(ubLHS[i],ubRHS[i]);
    for (int spt=0; spt<nSpts; spt++)
      for (int k=0; k<nFields; k++)
        eles[ic]->U_spts(spt,k) += unblankU(spt,k);
  }
#endif
}

void overComm::performProjection_static(vector<shared_ptr<ele>> &eles, vector<int> &eleMap, int order)
{
#ifndef _NO_MPI
  if (foundCells.size()<nproc)
    foundCells.resize(nproc);

  // Get the locations of the quadrature points for each target cell, and
  // the reference location of the points within the donor cells

  vector<matrix<double>> donorU(nproc);
  int nSpts = (order+1)*(order+1);
  if (nDims == 3) nSpts *= (order+1);
  int offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += foundCells[p-1].size();
    for (int i=0; i<foundCells[p].size(); i++) {
      for (int id=0; id<foundCellNDonors[p][i]; id++) {
        int ic = eleMap[foundCellDonors[p](i,id)];
        for (int spt=0; spt<nSpts; spt++)
          donorU[p].insertRow(&eles[ic]->U_spts(spt,0),INSERT_AT_END,nFields);
      }
    }
  }

  /* --- Perform the Galerkin Projection --- */

  // Perform integration for each target cell
  vector<matrix<double>> RHS(nproc);
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
      for (int ispt=0; ispt<nSpts; ispt++)
        for (int id=0; id<foundCellNDonors[p][i]; id++)
          for (int jspt=0; jspt<nSpts; jspt++)
            for (int k=0; k<nFields; k++)
              rhs(ispt,k) += massMatTDRow[(offset+i)*nSpts+ispt](id,jspt) * donorU[p](offsetD+id*nSpts+jspt,k);
      RHS[p].appendRows(rhs);
    }
  }

  // Send/recv the final target-cell data
  int strideR = nSpts*nFields;
  vector<matrix<double>> tmpUbRHS(nproc);
  sendRecvData(nCellsSend,nCellsRecv,RHS,tmpUbRHS,strideR);

  // Add contributions from superMeshes on each rank
  vector<matrix<double>> ubRHS(nUnblanks);
  for (auto &rhs:ubRHS) rhs.setup(nSpts,nFields);

  for (int p=0; p<nproc; p++) {
    if (p==rank) continue;
    for (int i=0; i<nCellsRecv[p]; i++) {
      for (int j=0; j<nSpts; j++) {
        for (int k=0; k<nFields; k++) {
          ubRHS[recvInds[p][i]](j,k) += tmpUbRHS[p](i,j*nFields+k);
        }
      }
    }
  }

  // Apply the new values to the unblank ele objects
  for (int i=0; i<nUnblanks; i++) {
    int ic = ubCells[i];
    for (int spt=0; spt<nSpts; spt++)
      for (int k=0; k<nFields; k++)
        eles[ic]->U_spts(spt,k) = 0;

    auto unblankU = solveCholesky(ubLHS[i],ubRHS[i]);
    for (int spt=0; spt<nSpts; spt++)
      for (int k=0; k<nFields; k++)
        eles[ic]->U_spts(spt,k) += unblankU(spt,k);
  }
#endif
}

vector<double> overComm::integrateErrOverset(vector<shared_ptr<ele>> &eles, map<int, oper> &opers, vector<int> &iblankCell, vector<int> &eleMap, int order, int quadOrder)
{
#ifndef _NO_MPI
  /* ---- Send Unblanked-Cell Nodes to All Grids ---- */

  int nOverlap = eles.size();

  int nv = (nDims==2) ? 4 : 8;

  vector<int> ubCells;
  Array<double,3> ubCellNodes(nOverlap,nv,nDims);
  int i = 0;

  for (int ie=0; ie<eles.size(); ie++) {
    if (iblankCell[eles[ie]->ID] == HOLE) continue;
    ubCells.push_back(ie);
    // Constraining this to just linear hexahedrons/quadrilaterals for the time being
    if (params->motion) {
      for (int j=0; j<nv; j++)
        for (int k=0; k<nDims; k++)
          ubCellNodes(i,j,k) = eles[ie]->nodesRK(j,k);
    } else {
      for (int j=0; j<nv; j++)
        for (int k=0; k<nDims; k++)
          ubCellNodes(i,j,k) = eles[ie]->nodes(j,k);
    }
    i++;
  }

  /* ---- Gather all cell bounding-box data on each grid ---- */

  int stride = nv*nDims;
  vector<int> nCells_rank;
  vector<double> ubNodes_rank;
  gatherData(nOverlap, stride, ubCellNodes.getData(), nCells_rank, ubNodes_rank);

  /* ---- Check Every Unblanked Cell for Donor Cells on This Grid ---- */

  // For use with ADT
  if (nDims == 2 and eleList.size() != eleMap.size()) {
    eleList.resize(eleMap.size());
    for (int i=0; i<eleMap.size(); i++) eleList[i] = i;
  }

  vector<vector<int>> foundCells(nproc); //.resize(nproc);
  vector<matrix<int>> foundCellDonors(nproc);
  vector<vector<int>> foundCellNDonors(nproc);
  vector<superMesh> supers(0);
  int offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += nCells_rank[p-1];

    if (gridIdList[p] == gridID) continue;

    foundCells[p].resize(0);
    foundCellDonors[p].setup(0,0);
    foundCellNDonors[p].resize(0);

    for (int i=0; i<nCells_rank[p]; i++) {
      // Get requested cell's bounding box [min & max extents]
      vector<double> targetBox(2*nDims);
      for (int dim=0; dim<nDims; dim++) {
        targetBox[dim]       =  1e15;
        targetBox[dim+nDims] = -1e15;
      }
      vector<point> targetNodes;
      for (int j=0; j<nv; j++) {
        point pt = point(&ubNodes_rank[(offset+i)*stride+nDims*j],nDims);
        targetNodes.push_back(pt);
        for (int dim=0; dim<nDims; dim++) {
          targetBox[dim]   = min(pt[dim],targetBox[dim]);
          targetBox[dim+nDims] = max(pt[dim],targetBox[dim+nDims]);
        }
      }

      // Find all possible donors using Tioga's ADT search (3D) or my brute-force search (2D)
      unordered_set<int> cellIDs;
      if (nDims == 2)
        adt->searchADT_box(eleList.data(),cellIDs,targetBox.data());
      else
        cellIDs = tg->findCellDonors(targetBox.data());

      // Remove interpolated field / fringe cells from integration
      set<int> fringeCells;
      for (auto &ic:cellIDs) {
        if (iblankCell[ic]!=NORMAL) {
          fringeCells.insert(ic);
        }
      }
      for (auto &ic:fringeCells) cellIDs.erase(ic);

      if (cellIDs.size() > 0) {
        vector<int> donorsIDs;

        foundCellNDonors[p].push_back(cellIDs.size());
        foundCells[p].push_back(i);
        for (auto &ic:cellIDs)
          donorsIDs.push_back(ic);

        // Setup the donor cells [on this grid] for the unblanked cell [on other grid]
        foundCellDonors[p].insertRowUnsized(donorsIDs);

        Array2D<point> donorPts;
        if (params->motion) {
          for (auto &ic:donorsIDs) {
            vector<point> tmpPts;
            for (uint npt = 0; npt < eles[eleMap[ic]]->nNodes; npt++)
              tmpPts.push_back(point(&eles[eleMap[ic]]->nodesRK(npt,0),nDims));
            donorPts.insertRow(tmpPts);
          }
        } else {
          for (auto &ic:donorsIDs) {
            vector<point> tmpPts;
            for (uint npt = 0; npt < eles[eleMap[ic]]->nNodes; npt++)
              tmpPts.push_back(point(&eles[eleMap[ic]]->nodes(npt,0),nDims));
            donorPts.insertRow(tmpPts);
          }
        }

        superMesh mesh(targetNodes,donorPts,quadOrder,nDims,rank,supers.size());
        supers.push_back(mesh);
      }
    }
  }

  /* --- Setup & Exchange Quadrature-Point Data --- */

  // Now that we have the local superMesh for each target, setup points for
  // use with numerical quadrature
  for (auto &mesh:supers)
    mesh.setupQuadrature();

  // Get the locations of the quadrature points for each target cell, and
  // interpolate the solution error to them
  vector<matrix<double>> superErr(supers.size());
  int nSpts = (order+1)*(order+1);
  if (nDims == 3) nSpts *= (order+1);
  offset = 0;
  for (int p=0; p<nproc; p++) {
    if (p>0) offset += foundCells[p-1].size();
    for (int i=0; i<foundCells[p].size(); i++) {
      vector<int> parents_tmp;
      matrix<double> qpts_tmp;
      supers[offset+i].getQpts(qpts_tmp,parents_tmp);

      for (int j=0; j<parents_tmp.size(); j++) {
        int ic = eleMap[foundCellDonors[p](i,parents_tmp[j])];
        point refLoc;
        bool isInEle = eles[ic]->getRefLocNelderMead(point(qpts_tmp[j],nDims),refLoc);

        if (!isInEle) {
          cout << "qpt: " << qpts_tmp(j,0) << ", " << qpts_tmp(j,1) << endl;
          auto box = eles[ic]->getBoundingBox();
          cout << "ele box: " << box[0] << ", " << box[1] << "; " << box[3] << ", " << box[4] << endl;
          cout << "ref loc: " << refLoc.x << ", " << refLoc.y << ", " << refLoc.z << endl;
          FatalError("Quadrature Point Reference Location not found in ele!");
        }

        vector<double> tmpU(nFields);
        int nSpts = eles[ic]->nSpts;
        matrix<double> tmpUspts(eles[ic]->nSpts,nFields);
        for (uint spt = 0; spt < nSpts; spt++)
          for (uint k = 0; k < nFields; k++)
            tmpUspts(spt,k) = eles[ic]->U_spts(spt,k);

        opers[eles[ic]->order].interpolateToPoint(tmpUspts, tmpU.data(), refLoc);
        vector<double> tmpErr = calcError(tmpU.data(),point(qpts_tmp[j],nDims),params);
        superErr[offset+i].insertRow(tmpErr);
      }
    }
  }

  /* --- Integrate the Entire Domain --- */

  vector<point> qpts;
  if (nDims == 2)
   qpts = getLocSpts(QUAD,quadOrder,string("Legendre"));
  else
   qpts = getLocSpts(HEX,quadOrder,string("Legendre"));

  auto wts = getQptWeights(quadOrder,nDims);

  matrix<double> quadPoints;
  for (auto &pt: qpts) quadPoints.insertRow({pt.x,pt.y,pt.z});

  vector<double> intErr(nFields);
  matrix<double> U_qpts;
  vector<double> detJac_qpts;
  for (uint ic=0; ic<eles.size(); ic++) {
    if (iblankCell[eles[ic]->ID]!=NORMAL) continue;

    int nSpts = eles[ic]->nSpts;
    matrix<double> tmpU_spts(eles[ic]->nSpts,nFields);
    vector<double> tmpDetJac_spts(eles[ic]->nSpts);
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint k = 0; k < nFields; k++)
        tmpU_spts(spt,k) = eles[ic]->U_spts(spt,k);
      tmpDetJac_spts[spt] = eles[ic]->detJac_spts(spt);
    }
    opers[eles[ic]->order].interpolateSptsToPoints(tmpU_spts, U_qpts, quadPoints);
    opers[eles[ic]->order].interpolateSptsToPoints(tmpDetJac_spts, detJac_qpts, quadPoints);
    for (uint i=0; i<qpts.size(); i++) {
      auto tmpErr = calcError(U_qpts[i], eles[ic]->calcPos(qpts[i]), params);
      for (int j=0; j<nFields; j++)
        intErr[j] += tmpErr[j] * wts[i] * detJac_qpts[i];
    }
  }

  /* --- Subtract the Overlap Region --- */

  for (int i=0; i<supers.size(); i++) {
    auto tmperr = supers[i].integrate(superErr[i]);

    if (tmperr.size()==0) continue; // Happens if supermesh is empty

    for (int j=0; j<nFields; j++)
      intErr[j] -= 0.5*tmperr[j];
  }

  vector<double> tmpErr = intErr;
  MPI_Allreduce(tmpErr.data(),intErr.data(),nFields,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (params->errorNorm == 2)
    for (auto &val:intErr) val = std::sqrt(std::abs(val));

  return intErr;
#endif
}

void overComm::exchangeOversetData(vector<shared_ptr<ele>> &eles, map<int, oper> &opers, vector<int> &eleMap)
{
#ifndef _NO_MPI
  U_out.resize(nproc);
  unordered_set<int> correctedEles;
  for (int p=0; p<nproc; p++) {
    if (params->oversetMethod == 1)
      U_out[p].setup(foundPts[p].size(),2*nFields);
    else
      U_out[p].setup(foundPts[p].size(),nFields);

    if (gridIdList[p] == gridID) continue;
    for (int i=0; i<foundPts[p].size(); i++) {
      point refPos = foundLocs[p][i];

      int ic = eleMap[foundEles[p][i]];

      uint nSpts = eles[ic]->nSpts;
      uint nFpts = eles[ic]->nFpts;

      if (ic<0 || ic>eles.size())
        FatalError("bad value of ic!");

      if (params->oversetMethod == 1) {
        /* Interpolate solution calculated from corrected flux
         * Need corrected flux for all donor cells, so calc discontinuous Fn & deltaFn
         *   for correction function scaling
         * Need to also keep track of whether flux is in ref/phys space and
         *   transform as required. */
        if (!correctedEles.count(ic)) {
          correctedEles.insert(ic);
          Array<double,3> tempF_spts(nDims,nSpts,nFields);
          matrix<double> norm_fpts(nFpts,nDims), tempFn_fpts(nFpts,nFields);
          vector<double> dA_fpts(nFpts);
          for (uint dim = 0; dim < nDims; dim++)
            for (uint spt = 0; spt < nSpts; spt++)
              for (uint k = 0; k < nFields; k++)
                tempF_spts(dim,spt,k) = eles[ic]->F_spts(dim,spt,k);
          for (uint fpt = 0; fpt < nFpts; fpt++) {
            dA_fpts[fpt] = eles[ic]->dA_fpts(fpt);
            for (uint dim = 0; dim < nDims; dim++)
              norm_fpts(fpt,dim) = eles[ic]->norm_fpts(fpt,dim);
          }
          if (params->motion) {
            opers[eles[ic]->order].applyExtrapolateFn(tempF_spts,norm_fpts,tempFn_fpts,dA_fpts);
          } else {
            opers[eles[ic]->order].applyExtrapolateFn(tempF_spts,tempFn_fpts);
          }

          for (uint fpt = 0; fpt < nFpts; fpt++)
            for (uint k = 0; k < nFields; k++)
              eles[ic]->disFn_fpts(fpt,k) = tempFn_fpts(fpt,k);
        }

        Array<double,3> tempF_spts;
        if (params->motion)
        {
          // Flux vector must be in ref. space in order to apply correction functions
          auto blah = eles[ic]->transformFlux_physToRef();
          for (uint dim = 0; dim < nDims; dim++)
            for (uint spt = 0; spt < nSpts; spt++)
              for (uint k = 0; k < nFields; k++)
                tempF_spts(dim,spt,k) = blah[dim](spt,k);
        }
        else
        {
          for (uint dim = 0; dim < nDims; dim++)
            for (uint spt = 0; spt < nSpts; spt++)
              for (uint k = 0; k < nFields; k++)
                tempF_spts(dim,spt,k) = eles[ic]->F_spts(dim,spt,k);
        }

        matrix<double> deltaFn(nFpts,nFields);
        for (uint fpt = 0; fpt < nFpts; fpt++)
          for (uint k = 0; k < nFields; k++)
            deltaFn(fpt,k) = eles[ic]->Fn_fpts(fpt,k) - eles[ic]->disFn_fpts(fpt,k);
        matrix<double> tempF_ref = opers[eles[ic]->order].interpolateCorrectedFlux(tempF_spts, deltaFn, refPos);

        vector<double> tempU(nFields);
        matrix<double> tempU_spts(nSpts,nFields);
        for (uint spt = 0; spt < nSpts; spt++)
          for (uint k = 0; k < nFields; k++)
            tempU_spts(spt,k) = eles[ic]->U_spts(spt,k);
        opers[eles[ic]->order].interpolateToPoint(tempU_spts, tempU.data(), refPos);

        // NOW we can transform flux vector back to physical space
        // [Recall: F_phys = (G/|G|) * F_ref]
        matrix<double> jacobian, invJaco;
        double detJac;
        eles[ic]->calcTransforms_point(jacobian,invJaco,detJac,refPos);

        matrix<double> tempF(nDims,nFields);
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            for (int k=0; k<nFields; k++)
              tempF(dim1,k) += jacobian(dim1,dim2) * tempF_ref(dim2,k) / detJac;

        if (params->motion) {
          for (int dim=0; dim<nDims; dim++)
            for (int k=0; k<nFields; k++)
              tempF(dim,k) += jacobian(dim,nDims) * tempU[k];
        }

        Vec3 outNorm = foundNorm[p][i];

        if (params->equation == NAVIER_STOKES) {
          vector<double> tempFn(nFields);
          for (int dim=0; dim<nDims; dim++)
            for (int field=0; field<nFields; field++)
              tempFn[field] += tempF(dim,field)*outNorm[dim];

          for (int k=0; k<nFields; k++) {
            U_out[p](i,k) = tempU[k];
            U_out[p](i,nFields+k) = tempFn[k];
          }
        }
        else if (params->equation == ADVECTION_DIFFUSION) {
          vector<double> tempFn(nFields);
          for (int dim=0; dim<nDims; dim++)
            for (int field=0; field<nFields; field++)
              tempFn[field] += tempF(dim,field)*outNorm[dim];

          // In case one of the advection speeds ~= 0, use max speed
          double vn = params->advectVx*outNorm.x + params->advectVy*outNorm.y;
          if (nDims == 3) vn += params->advectVz*outNorm.z;

          U_out[p](i,0) = tempFn[0]/vn;
        }
      }

      else {
        // Interpolate discontinuous solution [Original 'Artificial Boundary' Method]
        matrix<double> tempU_spts(nSpts,nFields);
        for (uint spt = 0; spt < nSpts; spt++)
          for (uint k = 0; k < nFields; k++)
            tempU_spts(spt,k) = eles[ic]->U_spts(spt,k);
        opers[eles[ic]->order].interpolateToPoint(tempU_spts, U_out[p][i], refPos);
      }

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

  // For flux-interp method, send both solution and normal flux
  int nVars = nFields;
  if (params->oversetMethod == 1) nVars *= 2;

  sendRecvData(nPtsSend,nPtsRecv,foundPts,recvPts,U_out,U_in,nVars,true);
#endif
}

void overComm::exchangeOversetGradient(vector<shared_ptr<ele>> &eles, map<int, oper> &opers, vector<int> &eleMap)
{
#ifndef _NO_MPI
  gradU_out.resize(nproc);

  for (int p=0; p<nproc; p++) {
    gradU_out[p].setup(foundPts[p].size(),nDims*nFields);
    if (gridIdList[p] == gridID) continue;

    for (int i=0; i<foundPts[p].size(); i++) {
      point refPos = foundLocs[p][i];
      int ic = eleMap[foundEles[p][i]];

      if (ic<0 || ic>eles.size())
        FatalError("bad value of ic!");

      // Interpolate corrected gradient
      // Need to also keep track of whether gradient is in ref/phys space and
      //   transform as required.
      vector<matrix<double>> tempDU_spts;
      uint nSpts = eles[ic]->nSpts;
      uint nFpts = eles[ic]->nFpts;
      if (params->motion) {
        // Gradient vector must be in ref. space in order to apply correction functions
        tempDU_spts = eles[ic]->transformGradU_physToRef();
      } else {
        for (uint dim = 0; dim < params->nDims; dim++)
          for (uint spt = 0; spt < eles[ic]->nSpts; spt++)
            for (uint k = 0; k < params->nFields; k++)
              tempDU_spts[dim](spt,k) = eles[ic]->dU_spts(dim,spt,k);
      }

      matrix<double> tempDU_ref(nDims,nFields);
      for (int dim=0; dim<nDims; dim++)
        opers[eles[ic]->order].interpolateToPoint(tempDU_spts[dim], tempDU_ref[dim], refPos);

      // NOW we can transform flux vector back to physical space
      // [Recall: F_phys = JGinv .dot. F_ref]
      matrix<double> jacobian, invJaco;
      double detJac;
      eles[ic]->calcTransforms_point(jacobian,invJaco,detJac,refPos);

      for (int dim1=0; dim1<nDims; dim1++)
        for (int dim2=0; dim2<nDims; dim2++)
          for (int k=0; k<nFields; k++)
            gradU_out[p](i,dim1*nFields+k) += invJaco(dim2,dim1) * tempDU_ref(dim2,k) / detJac;
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

  sendRecvData(nPtsSend,nPtsRecv,foundPts,recvPts,gradU_out,gradU_in,nDims*nFields,true);
#endif
}

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
