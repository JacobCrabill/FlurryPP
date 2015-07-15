/*!
 * \file output.cpp
 * \brief Functions for solution restart & visualization data output
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */
#include "../include/output.hpp"

#include <iomanip>

// Used for making sub-directories
#ifndef _NO_MPI
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
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
      e.updatePosSpts();
      e.updatePosFpts();
      e.setPpts();
    }
    for (uint spt=0; spt<e.getNSpts(); spt++) {
      V = e.getPrimitives(spt);
      pt = e.getPosSpt(spt);

      for (uint dim=0; dim<e.getNDims(); dim++) {
        dataFile << pt[dim] << ",";
      }
      if (e.getNDims() == 2) dataFile << "0.0,"; // output a 0 for z [2D]

      for (uint i=0; i<e.getNFields()-1; i++) {
        dataFile << V[i] << ",";
      }
      dataFile << V[e.getNFields()-1] << endl;
    }

    if (plotFpts) {
      for (uint fpt=0; fpt<e.getNFpts(); fpt++) {
        V = e.getPrimitivesFpt(fpt);
        pt = e.getPosFpt(fpt);

        for (uint dim=0; dim<e.getNDims(); dim++) {
          dataFile << pt[dim] << ",";
        }
        if (e.getNDims() == 2) dataFile << "0.0,"; // output a 0 for z [2D]

        for (uint i=0; i<e.getNFields()-1; i++) {
          dataFile << V[i] << ",";
        }
        dataFile << V[e.getNFields()-1] << endl;
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
  //if (params->meshType == OVERSET_MESH)
  //  sprintf(fileNameC,"%s_%.09d/%s%d_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],Solver->gridID,iter,Solver->gridRank);
  //else
    sprintf(fileNameC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],iter,params->rank);
#else
  /* --- Filename to write to --- */
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],iter);
#endif

  if (params->rank == 0)
    cout << "Writing ParaView file " << string(fileNameC) << "...  " << flush;

#ifndef _NO_MPI
  /* --- Write 'master' .pvtu file --- */
  if (params->rank == 0) {
    ofstream pVTU;
    char pvtuC[256];
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
      if (params->motion) {
        pVTU << "      <PDataArray type=\"Float32\" Name=\"GridVelocity\" NumberOfComponents=\"3\" />" << endl;
      }
    }
    pVTU << "    </PPointData>" << endl;
    pVTU << "    <PPoints>" << endl;
    pVTU << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" />" << endl;
    pVTU << "    </PPoints>" << endl;

    char filnameTmpC[256];
    for (int p=0; p<params->nproc; p++) {
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

  dataFile.open(fileNameC);
  dataFile.precision(16);

  // File header
  dataFile << "<?xml version=\"1.0\" ?>" << endl;
  dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
  dataFile << "	<UnstructuredGrid>" << endl;

  // If this is the initial file, need to extrapolate solution to flux points
  if (params->iter==params->initIter) Solver->extrapolateU();

  Solver->extrapolateUMpts();

  if (params->equation == NAVIER_STOKES) {
    if (params->calcEntropySensor)
      Solver->calcEntropyErr_spts();
    Solver->extrapolateSFpts();
    Solver->extrapolateSMpts();
  }

  for (auto& e:Solver->eles) {
    if (params->motion != 0) {
      e.updatePosSpts();
      e.updatePosFpts();
      e.setPpts();
    }

    // The combination of spts + fpts will be the plot points
    matrix<double> vPpts, gridVelPpts, errPpts;
    vector<point> ppts;
    e.getPrimitivesPlot(vPpts);
    if (params->motion)
      e.getGridVelPlot(gridVelPpts);
    ppts = e.getPpts();

    // Shock Capturing stuff
    double sensor;
    if(params->scFlag == 1) {
      sensor = e.getSensor();
    }

    if (params->equation == NAVIER_STOKES && params->calcEntropySensor)
      e.getEntropyErrPlot(errPpts);

    int nSubCells, nPpts;
    int nPpts1D = e.order+3;
    if (params->nDims == 2) {
      nSubCells = (e.order+2)*(e.order+2);
      nPpts = (e.order+3)*(e.order+3);
    }
    else if (params->nDims == 3) {
      nSubCells = (e.order+2)*(e.order+2)*(e.order+2);
      nPpts = (e.order+3)*(e.order+3)*(e.order+3);
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

    if(params->scFlag == 1){
      /* --- Shock Sensor --- */
      dataFile << "				<DataArray type=\"Float32\" Name=\"Sensor\" format=\"ascii\">" << endl;
      for(int k=0; k<nPpts; k++) {
        dataFile << sensor << " ";
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;
    }


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

      if (params->motion) {
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
    }

    if (params->equation == NAVIER_STOKES && params->calcEntropySensor) {
      /* --- Entropy Error Estimate --- */
      dataFile << "				<DataArray type=\"Float32\" Name=\"EntropyErr\" format=\"ascii\">" << endl;
      for(int k=0; k<nPpts; k++) {
        dataFile << std::abs(errPpts(k)) << " ";
      }
      dataFile << endl;
      dataFile << "				</DataArray>" << endl;
    }

    if (params->meshType == OVERSET_MESH) {
      /* --- TIOGA iBlank value --- */
      dataFile << "				<DataArray type=\"Float32\" Name=\"iBlank\" format=\"ascii\">" << endl;

      int iblank = 1;
      for (int i=0; i<Solver->Geo->c2nv[e.ID]; i++) {
        int iv = Solver->Geo->c2v(e.ID,i);
        iblank = min(iblank, Solver->Geo->iblank[iv]);
      }
      for(int k=0; k<nPpts; k++) {
        dataFile << iblank << " ";
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
#pragma omp parallel for
    for (uint e=0; e<Solver->eles.size(); e++) {
      auto resTmp = Solver->eles[e].getNormResidual(params->resType);
      if(checkNaN(resTmp)) FatalError("NaN Encountered in Solution Residual!");

      for (int i=0; i<params->nFields; i++)
        res[i] = max(res[i],resTmp[i]);
    }
  }
  else if (params->resType == 1 || params->resType == 2) {
    // 1-Norm or 2-Norm
#pragma omp parallel for
    for (uint e=0; e<Solver->eles.size(); e++) {
      auto resTmp = Solver->eles[e].getNormResidual(params->resType);
      if(checkNaN(resTmp)) FatalError("NaN Encountered in Solution Residual!");

      for (int i=0; i<params->nFields; i++)
        res[i] += resTmp[i];
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

  // If taking 2-norm, res is sum squared; take sqrt to complete
  if (params->rank == 0) {
    if (params->resType == 2) {
      for (auto& R:res) R = sqrt(R);
    }

    int colW = 16;
    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    if (iter==1 || (iter/params->monitorResFreq)%25==0) {
      cout << endl;
      cout << setw(8) << left << "Iter";
      if (params->equation == ADVECTION_DIFFUSION) {
        cout << " Residual " << endl;
      }else if (params->equation == NAVIER_STOKES) {
        cout << setw(colW) << left << "rho";
        cout << setw(colW) << left << "rhoU";
        cout << setw(colW) << left << "rhoV";
        if (params->nDims == 3)
          cout << setw(colW) << left << "rhoW";
        cout << setw(colW) << left << "rhoE";
        if (params->dtType == 1)
          cout << setw(colW) << left << "deltaT";
      }
      cout << endl;
    }

    cout << setw(8) << left << iter;
    for (int i=0; i<params->nFields; i++) {
      cout << setw(colW) << left << res[i];
    }
    if (params->dtType == 1)
      cout << setw(colW) << left << params->dt;
    cout << endl;
  }
}

void writeMeshTecplot(solver* Solver, input* params)
{
  if (params->nDims == 2) return;

  ofstream dataFile;

  char fileNameC[100];
  string fileName = params->dataFileName;

  geo* Geo = Solver->Geo;

#ifndef _NO_MPI
  /* --- All processors write their solution to their own .vtu file --- */
  sprintf(fileNameC,"%s/%s_%d.plt",&fileName[0],&fileName[0],params->rank);
#else
  /* --- Filename to write to --- */
  sprintf(fileNameC,"%s.plt",&fileName[0]);
#endif

  //if (params->rank == 0)
  //  cout << "Writing Tecplot mesh file " << string(fileNameC) << "...  " << flush;
  cout << "Writing Tecplot mesh file " << string(fileNameC) << "...  " << endl;

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

  cout << "rank " << params->rank << ": Opening file for writing" << endl;

  dataFile.open(fileNameC);
  dataFile.precision(16);

  // Count wall-boundary nodes
  int nNodesWall = Geo->iwall.size();
  // Count overset-boundary nodes
  int nNodesOver = Geo->iover.size();

  int gridID = Geo->gridID;

  cout << "nNodesWall = " << nNodesWall << ", nNodesOver = " << nNodesOver << ", gridID = " << gridID << endl;

  int nPrism = 0;
  int nNodes = Geo->nVerts;
  int nCells = Geo->nEles;
  int nHex = nCells;

  cout << "nNodes = " << nNodes << ", nEles = " << nCells << endl;

  dataFile << "# " << nPrism << " " << nHex << " " << nNodes << " " << nCells << " " << nNodesWall << " " << nNodesOver << endl;
  dataFile << "TITLE = \"" << fileName << "\"" << endl;
  dataFile << "VARIABLES = \"X\", \"Y\", \"Z\", \"bodyTag\", \"IBLANK\"" << endl;
  dataFile << "ZONE T = \"VOL_MIXED\", N=" << nNodes << ", E=" << nCells << ", ET=BRICK, F=FEPOINT" << endl;

  for (int iv=0; iv<nNodes; iv++) {
    //cout << Geo->xv(iv,0) << " " << Geo->xv(iv,1) << " " << Geo->xv(iv,2) << " " << gridID << endl;
    dataFile << Geo->xv(iv,0) << " " << Geo->xv(iv,1) << " " << Geo->xv(iv,2) << " " << gridID << " " << Geo->iblank[iv] << endl;
  }

  for (int ic=0; ic<nCells; ic++) {
    for (int j=0; j<8; j++) {
      dataFile << Geo->c2v(ic,j)+1 << " ";
    }
    dataFile << endl;
  }

  // output wall-boundary node IDs
  for (auto& iv:Geo->iwall)
    dataFile << iv+1 << endl;

  // output overset-boundary node IDs
  for (auto& iv:Geo->iover)
    dataFile << iv+1 << endl;

  dataFile.close();

  if (params->rank == 0) cout << "done." << endl;
}
