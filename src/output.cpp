/*!
 * \file output.cpp
 * \brief Functions for solution restart & visualization data output
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
#include "output.hpp"

#include <iomanip>
#include <string>

// Used for making sub-directories (for MPI and 'time-stamp' files)
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef _NO_MPI
#include "mpi.h"
#endif

void writeData(solver *Solver, input *params)
{
  if (params->plotType == 0) {
    writeCSV(Solver,params);
  }
  else if (params->plotType == 1) {
    writeParaview(Solver,params);
  }

  /* Write out mesh in Tecplot format, with IBLANK data [Overset cases only] */
  if (params->meshType==OVERSET_MESH && params->writeIBLANK)
    writeMeshTecplot(Solver,params);
}

void writeCSV(solver *Solver, input *params)
{
  ofstream dataFile;
  int iter = params->iter;

  char fileNameC[256];
  string fileName = params->dataFileName;
  sprintf(fileNameC,"%s.csv.%.09d",&fileName[0],iter);

  dataFile.precision(15);
  dataFile.setf(ios_base::fixed);

  dataFile.open(fileNameC);

  // Vector of primitive variables
  vector<double> V;
  // Location of solution point
  point pt;

  // Write header:
  // x  y  z(=0)  rho  [u  v  p]
  dataFile << "x,y,z,";
  if (params->equation == ADVECTION_DIFFUSION) {
    dataFile << "rho" << endl;
  }
  else if (params->equation == NAVIER_STOKES) {
    dataFile << "rho,u,v,p" << endl;
  }

  bool plotFpts = true;

  if (plotFpts)
    Solver->extrapolateU();

  // Solution data
  for (auto& e:Solver->eles) {
    if (params->motion != 0) {
      e->updatePosSpts();
      e->updatePosFpts();
      e->setPpts();
    }
    for (uint spt=0; spt<e->getNSpts(); spt++) {
      V = e->getPrimitives(spt);
      pt = e->getPosSpt(spt);

      for (uint dim=0; dim<e->getNDims(); dim++) {
        dataFile << pt[dim] << ",";
      }
      if (e->getNDims() == 2) dataFile << "0.0,"; // output a 0 for z [2D]

      for (uint i=0; i<e->getNFields()-1; i++) {
        dataFile << V[i] << ",";
      }
      dataFile << V[e->getNFields()-1] << endl;
    }

    if (plotFpts) {
      for (uint fpt=0; fpt<e->getNFpts(); fpt++) {
        V = e->getPrimitivesFpt(fpt);
        pt = e->getPosFpt(fpt);

        for (uint dim=0; dim<e->getNDims(); dim++) {
          dataFile << pt[dim] << ",";
        }
        if (e->getNDims() == 2) dataFile << "0.0,"; // output a 0 for z [2D]

        for (uint i=0; i<e->getNFields()-1; i++) {
          dataFile << V[i] << ",";
        }
        dataFile << V[e->getNFields()-1] << endl;
      }
    }
  }

  dataFile.close();
}

void writeParaview(solver *Solver, input *params)
{
  ofstream dataFile;
  int iter = params->iter;

  char fileNameC[256];
  string fileName = params->dataFileName;

#ifndef _NO_MPI
  /* --- All processors write their solution to their own .vtu file --- */
  if (params->meshType == OVERSET_MESH)
    sprintf(fileNameC,"%s_%.09d/%s%d_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],Solver->gridID,iter,Solver->gridRank);
  else
    sprintf(fileNameC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],iter,params->rank);
#else
  /* --- Filename to write to --- */
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],iter);
#endif

  char Iter[10];
  sprintf(Iter,"%.09d",iter);

  if (params->rank == 0)
    cout << "Writing ParaView file " << params->dataFileName << "_" << string(Iter) << ".vtu...  " << flush;

#ifndef _NO_MPI
  /* --- Write 'master' .pvtu file (for each grid, if overset) --- */
  if (Solver->gridRank == 0) {
    ofstream pVTU;
    char pvtuC[256];
    if (params->meshType == OVERSET_MESH)
      sprintf(pvtuC,"%s%d_%.09d.pvtu",&fileName[0],Solver->gridID,iter);
    else
      sprintf(pvtuC,"%s_%.09d.pvtu",&fileName[0],iter);

    pVTU.open(pvtuC);

    pVTU << "<?xml version=\"1.0\" ?>" << endl;
    pVTU << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
    pVTU << "  <PUnstructuredGrid GhostLevel=\"1\">" << endl;
    // NOTE: Must be careful with order here [particularly of vector data], or else ParaView gets confused
    pVTU << "    <PPointData Scalars=\"Density\" Vectors=\"Velocity\" >" << endl;
    pVTU << "      <PDataArray type=\"Float32\" Name=\"Density\" />" << endl;
    if (params->equation == NAVIER_STOKES) {
      pVTU << "      <PDataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" />" << endl;
      pVTU << "      <PDataArray type=\"Float32\" Name=\"Pressure\" />" << endl;
      if (params->calcEntropySensor) {
        pVTU << "      <PDataArray type=\"Float32\" Name=\"EntropyErr\" />" << endl;
      }
    }
    if (params->scFlag == 1) {
      pVTU << "      <PDataArray type=\"Float32\" Name=\"Sensor\" />" << endl;
    }
    if (params->motion) {
      pVTU << "      <PDataArray type=\"Float32\" Name=\"GridVelocity\" NumberOfComponents=\"3\" />" << endl;
    }
    if (params->meshType == OVERSET_MESH && params->writeIBLANK) {
      pVTU << "      <PDataArray type=\"Float32\" Name=\"IBLANK\" />" << endl;
    }
    pVTU << "    </PPointData>" << endl;
    pVTU << "    <PPoints>" << endl;
    pVTU << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" />" << endl;
    pVTU << "    </PPoints>" << endl;

    char filnameTmpC[256];
    for (int p=0; p<Solver->nprocPerGrid; p++) {
      if (params->meshType == OVERSET_MESH)
        sprintf(filnameTmpC,"%s_%.09d/%s%d_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],Solver->gridID,iter,p);
      else
        sprintf(filnameTmpC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],iter,p);
      pVTU << "    <Piece Source=\"" << string(filnameTmpC) << "\" />" << endl;
    }
    pVTU << "  </PUnstructuredGrid>" << endl;
    pVTU << "</VTKFile>" << endl;

    pVTU.close();

    char datadirC[256];
    char *datadir = &datadirC[0];
    sprintf(datadirC,"%s_%.09d",&fileName[0],iter);

    /* --- Master node creates a subdirectory to store .vtu files --- */
    if (params->rank == 0) {
      struct stat st = {0};
      if (stat(datadir, &st) == -1) {
        mkdir(datadir, 0755);
      }
    }
  }

  /* --- Wait for all processes to get here, otherwise there won't be a
   *     directory to put .vtus into --- */
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* --- There's no possible way to put the simulation time into a ParaView file,
         so we have to create a separate file to store the time values --- */
  char datadirC[256];
  char *datadir = &datadirC[0];
  sprintf(datadirC,"%s_time",&fileName[0]);
  char timeFileC[256];
  sprintf(timeFileC,"%s_time/%d",&fileName[0],iter);

  if (params->rank == 0) {
    struct stat st = {0};
    if (stat(datadir, &st) == -1) {
      mkdir(datadir, 0755);
    }
    dataFile.open(timeFileC);
    dataFile << params->time << endl;
    dataFile.clear();
    dataFile.close();
  }

  // Write the cell iblank data for restarting purposes
  if (params->meshType == OVERSET_MESH) {
    for (int p=0; p<params->nproc; p++) {
      if (p == params->rank) {
        dataFile.open(timeFileC,ios::app);
        dataFile << params->rank << " ";
        for (int i=0; i<Solver->Geo->nEles; i++) {
          if (Solver->Geo->eleMap[i]<0)
            dataFile << HOLE << " ";
          else
            dataFile << Solver->Geo->iblankCell[i] << " ";
        }
        dataFile << endl;
        dataFile.clear();
        dataFile.close();
      }
#ifndef _NO_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  }

  dataFile.open(fileNameC);
  dataFile.precision(16);

  // File header
  dataFile << "<?xml version=\"1.0\" ?>" << endl;
  dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
  dataFile << "	<UnstructuredGrid>" << endl;

  // If this is the initial file, need to extrapolate solution to flux points
  if (params->iter==params->initIter) Solver->extrapolateU();

  Solver->extrapolateUMpts();

  if (params->motion && params->nDims==3)
    Solver->extrapolateGridVelMpts();

  if (params->equation == NAVIER_STOKES) {
    if (params->squeeze) {
      Solver->calcAvgSolution();
      Solver->checkEntropyPlot();
    }

    if (params->calcEntropySensor) {
      Solver->calcEntropyErr_spts();
      Solver->extrapolateSFpts();
      Solver->extrapolateSMpts();
    }
  }

  for (auto& e:Solver->eles) {
    if (params->meshType == OVERSET_MESH && Solver->Geo->iblankCell[e->ID]!=NORMAL) continue;

    if (params->motion != 0) {
      e->updatePosSpts();
      e->updatePosFpts();
      e->setPpts();
    }

    // The combination of spts + fpts will be the plot points
    matrix<double> vPpts, gridVelPpts, errPpts;
    vector<point> ppts;
    e->getPrimitivesPlot(vPpts);
    if (params->motion)
      e->getGridVelPlot(gridVelPpts);
    ppts = e->getPpts();

    // Shock Capturing stuff
    double sensor;
    if(params->scFlag == 1) {
      sensor = e->getSensor();
    }

    if (params->equation == NAVIER_STOKES && params->calcEntropySensor)
      e->getEntropyErrPlot(errPpts);

    int nSubCells, nPpts;
    int nPpts1D = e->order+3;
    if (params->nDims == 2) {
      nSubCells = (e->order+2)*(e->order+2);
      nPpts = (e->order+3)*(e->order+3);
    }
    else if (params->nDims == 3) {
      nSubCells = (e->order+2)*(e->order+2)*(e->order+2);
      nPpts = (e->order+3)*(e->order+3)*(e->order+3);
    }
    else
      FatalError("Invalid dimensionality [nDims].");

    // Write cell header
    dataFile << "		<Piece NumberOfPoints=\"" << nPpts << "\" NumberOfCells=\"" << nSubCells << "\">" << endl;

    /* ==== Write out solution to file ==== */

    dataFile << "			<PointData>" << endl;

    /* --- Density --- */
    dataFile << "				<DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">" << endl;
    for(int k=0; k<nPpts; k++) {
      dataFile << vPpts(k,0) << " ";
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    if (params->equation == NAVIER_STOKES) {
      /* --- Velocity --- */
      dataFile << "				<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">" << endl;
      for(int k=0; k<nPpts; k++) {
        dataFile << vPpts(k,1) << " " << vPpts(k,2) << " ";

        // In 2D the z-component of velocity is not stored, but Paraview needs it so write a 0.
        if(params->nDims==2) {
          dataFile << 0.0 << " ";
        }
        else {
          dataFile << vPpts(k,3) << " ";
        }
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;

      /* --- Pressure --- */
      dataFile << "				<DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">" << endl;
      for(int k=0; k<nPpts; k++) {
        dataFile << vPpts(k,params->nDims+1) << " ";
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;

      if (params->calcEntropySensor) {
        /* --- Entropy Error Estimate --- */
        dataFile << "				<DataArray type=\"Float32\" Name=\"EntropyErr\" format=\"ascii\">" << endl;
        for(int k=0; k<nPpts; k++) {
          dataFile << std::abs(errPpts(k)) << " ";
        }
        dataFile << endl;
        dataFile << "				</DataArray>" << endl;
      }
    }

    if(params->scFlag == 1) {
      /* --- Shock Sensor --- */
      dataFile << "				<DataArray type=\"Float32\" Name=\"Sensor\" format=\"ascii\">" << endl;
      for(int k=0; k<nPpts; k++) {
        dataFile << sensor << " ";
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;
    }

    if (params->motion > 0) {
      /* --- Grid Velocity --- */
      dataFile << "				<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"GridVelocity\" format=\"ascii\">" << endl;
      for(int k=0; k<nPpts; k++) {
        // Divide momentum components by density to obtain velocity components
        dataFile << gridVelPpts(k,0) << " " << gridVelPpts(k,1) << " ";

        // In 2D the z-component of velocity is not stored, but Paraview needs it so write a 0.
        if(params->nDims==2) {
          dataFile << 0.0 << " ";
        }
        else {
          dataFile << gridVelPpts(k,2) << " ";
        }
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;
    }

    if (params->meshType == OVERSET_MESH && params->writeIBLANK) {
      /* --- TIOGA iBlank value --- */
      dataFile << "				<DataArray type=\"Float32\" Name=\"IBLANK\" format=\"ascii\">" << endl;

      for(int k=0; k<nPpts; k++) {
        dataFile << Solver->Geo->iblankCell[e->ID] << " ";
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;
    }

    /* --- End of Cell's Solution Data --- */

    dataFile << "			</PointData>" << endl;

    /* ==== Write Out Cell Points & Connectivity==== */

    /* --- Write out the plot point coordinates --- */
    dataFile << "			<Points>" << endl;
    dataFile << "				<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

    // Loop over plot points in element
    for(int k=0; k<nPpts; k++) {
      for(int l=0;l<params->nDims;l++) {
        dataFile << ppts[k][l] << " ";
      }

      // If 2D, write a 0 as the z-component
      if(params->nDims == 2) {
        dataFile << "0 ";
      }
    }

    dataFile << endl;
    dataFile << "				</DataArray>" << endl;
    dataFile << "			</Points>" << endl;

    /* --- Write out Cell data: connectivity, offsets, element types --- */
    dataFile << "			<Cells>" << endl;

    /* --- Write connectivity array --- */
    dataFile << "				<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;

    if (params->nDims == 2) {
      for (int j=0; j<nPpts1D-1; j++) {
        for (int i=0; i<nPpts1D-1; i++) {
          dataFile << j*nPpts1D     + i   << " ";
          dataFile << j*nPpts1D     + i+1 << " ";
          dataFile << (j+1)*nPpts1D + i+1 << " ";
          dataFile << (j+1)*nPpts1D + i   << " ";
          dataFile << endl;
        }
      }
    }
    else if (params->nDims == 3) {
      for (int k=0; k<nPpts1D-1; k++) {
        for (int j=0; j<nPpts1D-1; j++) {
          for (int i=0; i<nPpts1D-1; i++) {
            dataFile << i   + nPpts1D*(j   + nPpts1D*k) << " ";
            dataFile << i+1 + nPpts1D*(j   + nPpts1D*k) << " ";
            dataFile << i+1 + nPpts1D*(j+1 + nPpts1D*k) << " ";
            dataFile << i   + nPpts1D*(j+1 + nPpts1D*k) << " ";

            dataFile << i   + nPpts1D*(j   + nPpts1D*(k+1)) << " ";
            dataFile << i+1 + nPpts1D*(j   + nPpts1D*(k+1)) << " ";
            dataFile << i+1 + nPpts1D*(j+1 + nPpts1D*(k+1)) << " ";
            dataFile << i   + nPpts1D*(j+1 + nPpts1D*(k+1)) << " ";

            dataFile << endl;
          }
        }
      }
    }
    dataFile << "				</DataArray>" << endl;

    // Write cell-node offsets
    int nvPerCell;
    if (params->nDims == 2) nvPerCell = 4;
    else                    nvPerCell = 8;
    dataFile << "				<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    for(int k=0; k<nSubCells; k++){
      dataFile << (k+1)*nvPerCell << " ";
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    // Write VTK element type
    // 5 = tri, 9 = quad, 10 = tet, 12 = hex
    int eType;
    if (params->nDims == 2) eType = 9;
    else                    eType = 12;
    dataFile << "				<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    for(int k=0; k<nSubCells; k++) {
      dataFile << eType << " ";
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    /* --- Write cell and piece footers --- */
    dataFile << "			</Cells>" << endl;
    dataFile << "		</Piece>" << endl;
  }

  /* --- Write footer of file & close --- */
  dataFile << "	</UnstructuredGrid>" << endl;
  dataFile << "</VTKFile>" << endl;

  dataFile.close();

  if (params->rank == 0) cout << "done." <<  endl;
}


void writeResidual(solver *Solver, input *params)
{
  vector<double> res(params->nFields);
  int iter = params->iter;

  if (params->resType == 3) {
    // Infinity Norm
    for (uint e=0; e<Solver->eles.size(); e++) {
      if (params->meshType == OVERSET_MESH && Solver->Geo->iblankCell[Solver->eles[e]->ID]!=NORMAL) continue;
      auto resTmp = Solver->eles[e]->getNormResidual(params->resType);
      if(checkNaN(resTmp)) {
        cout << "rank " << params->rank << ", ele " << e << ": ";
        auto box = Solver->eles[e]->getBoundingBox();
        cout << " minPt = " << box[0] << "," << box[1] << "," << box[2] << ", maxPt = " << box[3] << "," << box[4] << "," << box[5] << endl;
        FatalError("NaN Encountered in Solution Residual!");
      }

      for (int i=0; i<params->nFields; i++)
        res[i] = max(res[i],resTmp[i]);
    }
  }
  else if (params->resType == 1 || params->resType == 2) {
    // 1-Norm or 2-Norm
    for (uint e=0; e<Solver->eles.size(); e++) {
      if (params->meshType == OVERSET_MESH && Solver->Geo->iblankCell[Solver->eles[e]->ID]!=NORMAL) continue;
      auto resTmp = Solver->eles[e]->getNormResidual(params->resType);
      if(checkNaN(resTmp)) {
        cout << "rank " << params->rank << ", ele " << e << ": " << flush;
        auto box = Solver->eles[e]->getBoundingBox();
        cout << " minPt = " << box[0] << "," << box[1] << "," << box[2] << ", maxPt = " << box[3] << "," << box[4] << "," << box[5] << endl;
        FatalError("NaN Encountered in Solution Residual!");
      }

      for (int i=0; i<params->nFields; i++)
        res[i] += resTmp[i];
    }
  }

  vector<double> force(6);
  if (params->equation == NAVIER_STOKES) {
    auto fTmp = Solver->computeWallForce();
#ifndef _NO_MPI
    MPI_Reduce(fTmp.data(), force.data(), 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    force = fTmp;
#endif
    for (auto &f:force) f /= (0.5*params->rhoBound*params->Uinf*params->Uinf);

    double alpha = std::atan2(params->vBound,params->uBound);
    fTmp = force;
    force[0] = fTmp[0]*cos(alpha) + fTmp[1]*sin(alpha);  // Rotate to align with freestream
    force[1] = fTmp[1]*cos(alpha) - fTmp[0]*sin(alpha);
    if (params->viscous) {
      force[3] = fTmp[3]*cos(alpha) + fTmp[4]*sin(alpha);
      force[4] = fTmp[4]*cos(alpha) - fTmp[3]*sin(alpha);
    }
  }

#ifndef _NO_MPI
  if (params->nproc > 1) {
    if (params->resType == 3) {
      vector<double> resTmp = res;
      MPI_Reduce(resTmp.data(), res.data(), params->nFields, MPI_DOUBLE, MPI_MAX, 0,MPI_COMM_WORLD);
    }
    else if (params->resType == 1 || params->resType == 2) {
      vector<double> resTmp = res;
      MPI_Reduce(resTmp.data(), res.data(), params->nFields, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    }
  }
#endif

  if (params->rank == 0) {
    // If taking 2-norm, res is sum squared; take sqrt to complete
    if (params->resType == 2) {
      for (auto& R:res) R = sqrt(R);
    }

    /* --- Print the residual and force coefficients in the terminal --- */

    int colW = 16;
    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    if (iter==params->initIter+1 || (iter/params->monitorResFreq)%25==0) {
      cout << endl;
      cout << setw(8) << left << "Iter" << "Var  ";
      if (params->equation == ADVECTION_DIFFUSION) {
        cout << "Residual" << endl;
      }else if (params->equation == NAVIER_STOKES) {
        cout << setw(colW) << left << "rho";
        cout << setw(colW) << left << "rhoU";
        cout << setw(colW) << left << "rhoV";
        if (params->nDims == 3)
          cout << setw(colW) << left << "rhoW";
        cout << setw(colW) << left << "rhoE";
        if (params->dtType == 1)
          cout << setw(colW) << left << "deltaT";
        cout << setw(colW) << left << "CD";
        cout << setw(colW) << left << "CL";
        if (params->nDims == 3)
          cout << setw(colW) << left << "CN";
      }
      cout << endl;
    }

    // Print residuals
    cout << setw(8) << left << iter << "Res  ";
    for (int i=0; i<params->nFields; i++) {
      cout << setw(colW) << left << res[i];
    }

    // Print time step (for CFL time-stepping)
    if (params->dtType == 1)
      cout << setw(colW) << left << params->dt;

    // Print wall force coefficients
    if (params->equation == NAVIER_STOKES) {
      for (int dim=0; dim<params->nDims; dim++) {
        cout << setw(colW) << left << force[dim]+force[3+dim];
      }
    }

    cout << endl;

    /* --- Write the residual and force coefficients to the history file --- */

    ofstream histFile;
    string fileName = params->dataFileName + ".hist";
    histFile.open(fileName.c_str(),ofstream::app);

    histFile.precision(5);
    histFile.setf(ios::scientific, ios::floatfield);
    if (iter==params->initIter+1 || (iter/params->monitorResFreq)%25==0) {
      histFile << endl;
      histFile << setw(8) << left << "Iter";
      histFile << setw(colW) << left << "Time";
      histFile << "Var  ";
      if (params->equation == ADVECTION_DIFFUSION) {
        histFile << setw(colW) << left << "Residual" << endl;
      }else if (params->equation == NAVIER_STOKES) {
        histFile << setw(colW) << left << "rho";
        histFile << setw(colW) << left << "rhoU";
        histFile << setw(colW) << left << "rhoV";
        if (params->nDims == 3)
          histFile << setw(colW) << left << "rhoW";
        histFile << setw(colW) << left << "rhoE";
        if (params->dtType == 1)
          histFile << setw(colW) << left << "deltaT";
        histFile << setw(colW) << left << "CDinv";
        histFile << setw(colW) << left << "CLinv";
        if (params->nDims == 3)
          histFile << setw(colW) << left << "CNinv";
        if (params->viscous) {
          histFile << setw(colW) << left << "CDvis";
          histFile << setw(colW) << left << "CLvis";
          if (params->nDims == 3)
            histFile << setw(colW) << left << "CNvis";
        }
        histFile << setw(colW) << left << "CDtot";
        histFile << setw(colW) << left << "CLtot";
        if (params->nDims == 3)
          histFile << setw(colW) << left << "CNtot";
      }
      histFile << endl;
    }

    // Write residuals
    histFile << setw(8) << left << iter;
    histFile << setw(colW) << left << params->time;
    histFile << "Res  ";
    for (int i=0; i<params->nFields; i++) {
      histFile << setw(colW) << left << res[i];
    }

    // Write time step (for CFL time-stepping)
    if (params->dtType == 1)
      histFile << setw(colW) << left << params->dt;

    // Write inviscid wall force coefficients
    if (params->equation == NAVIER_STOKES) {
      for (int dim=0; dim<params->nDims; dim++)
        histFile << setw(colW) << left << force[dim];              // Convective force coeffs.
      if (params->viscous)
        for (int dim=0; dim<params->nDims; dim++)
          histFile << setw(colW) << left << force[3+dim];          // Viscous force coeffs.
      for (int dim=0; dim<params->nDims; dim++)
        histFile << setw(colW) << left << force[dim]+force[3+dim]; // Total force coeffs.
    }

    histFile << endl;
    histFile.close();
  }
}

void writeAllError(solver *Solver, input *params)
{
  if (params->testCase > 0) {
    params->errorNorm = 0;
    if (params->rank == 0)
      cout << "Integrated conservation error:" << endl;
    writeError(Solver,params);

    params->errorNorm = 1;
    if (params->rank == 0)
      cout << "Integral L1 error:" << endl;
    writeError(Solver,params);

    params->errorNorm = 2;
    if (params->rank == 0)
      cout << "Integral L2 error:" << endl;
    writeError(Solver,params);
  } else {
    params->errorNorm = 0;
    if (params->rank == 0)
      cout << "Integrated conservative variables:" << endl;
    writeError(Solver,params);
  }
}

void writeError(solver *Solver, input *params)
{
  if (params->meshType != OVERSET_MESH && !params->testCase) return;

  // For implemented test cases, calculcate the L1 error over the overset domain

  vector<double> err;

  if (params->meshType == OVERSET_MESH)
    err = Solver->integrateErrorOverset();
  else
    err = Solver->integrateError();

  if (params->rank == 0) {
    /* --- Write the error out to the terminal --- */

    int colw = 16;
    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);

    cout << setw(8) << left << params->iter << "Err  ";
    for (int i=0; i<err.size(); i++)
      cout << setw(colw) << left << std::abs(err[i]);
    cout << endl;

    /* --- Write the error out to the history file --- */

    ofstream histFile;
    string fileName = params->dataFileName + ".hist";
    histFile.open(fileName.c_str(),ofstream::app);

    histFile.precision(5);
    histFile.setf(ios::scientific, ios::floatfield);

    histFile << setw(8) << left << params->iter;
    histFile << setw(colw) << left << params->time;
    histFile << "Err  ";
    for (int i=0; i<params->nFields; i++) {
      histFile << setw(colw) << left << std::abs(err[i]);
    }
    histFile << endl;
    histFile.close();
  }
}

void writeMeshTecplot(solver* Solver, input* params)
{
  ofstream dataFile;

  char fileNameC[100];
  string fileName = params->dataFileName;

  geo* Geo = Solver->Geo;

#ifndef _NO_MPI
  /* --- All processors write their solution to their own .vtu file --- */
  sprintf(fileNameC,"%s/%s_%d_%d.plt",&fileName[0],&fileName[0],params->iter,params->rank);
#else
  /* --- Filename to write to --- */
  sprintf(fileNameC,"%s.plt",&fileName[0]);
#endif

  if (params->rank == 0)
    cout << "Writing Tecplot mesh file " << string(fileNameC) << "...  " << flush;

#ifndef _NO_MPI
  /* --- Write folder for output mesh files --- */
  if (params->rank == 0) {

    char datadirC[100];
    char *datadir = &datadirC[0];
    sprintf(datadirC,"%s",&fileName[0]);

    /* --- Master node creates a subdirectory to store .vtu files --- */
    if (params->rank == 0) {
      struct stat st = {0};
      if (stat(datadir, &st) == -1) {
        mkdir(datadir, 0755);
      }
    }
  }

  /* --- Wait for all processes to get here, otherwise there won't be a
   *     directory to put .vtus into --- */
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  dataFile.open(fileNameC);
  dataFile.precision(16);

  // Count wall-boundary nodes
  int nNodesWall = Geo->iwall.size();
  // Count overset-boundary nodes
  int nNodesOver = Geo->iover.size();

  int gridID = Geo->gridID;

  int nPrism = 0;
  int nNodes = Geo->nVerts;
  int nCells = Geo->nEles;
  int nHex = nCells;

  dataFile << "# " << nPrism << " " << nHex << " " << nNodes << " " << nCells << " " << nNodesWall << " " << nNodesOver << endl;
  dataFile << "TITLE = \"" << fileName << "\"" << endl;
  dataFile << "VARIABLES = \"X\", \"Y\", \"Z\", \"bodyTag\", \"IBLANK\", \"IBLANKCELL\"" << endl;
  string ET;
  if (params->nDims==2)
    ET = "QUADRILATERAL";
  else
    ET = "BRICK";
  dataFile << "ZONE T = \"VOL_MIXED\", N=" << nCells*4 << ", E=" << nCells << ", ET=" << ET << ", F=FEPOINT" << endl;

  for (int ic=0; ic<nCells; ic++) {
    for (int j=0; j<4; j++) {
      dataFile << Geo->xv(Geo->c2v(ic,j),0) << " " << Geo->xv(Geo->c2v(ic,j),1) << " " << 0.0 << " " << gridID << " " << Geo->iblank[Geo->c2v(ic,j)] << " " << Geo->iblankCell[ic] << endl;
    }
  }

  int nv = (params->nDims==2) ? 4 : 8; // Ignoring edge nodes for quadratic elements
  for (int ic=0; ic<nCells; ic++) {
    for (int j=0; j<nv; j++) {
      dataFile << ic*nv + j + 1 << " ";
    }
    dataFile << endl;
  }

//  // output wall-boundary node IDs
//  for (auto& iv:Geo->iwall)
//    dataFile << iv+1 << endl;

//  // output overset-boundary node IDs
//  for (auto& iv:Geo->iover)
//    dataFile << iv+1 << endl;

  dataFile.close();

  if (params->rank == 0) cout << "done." << endl;
}
