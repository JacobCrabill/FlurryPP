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
  if (params->plot_type == 0) {
    writeCSV(Solver,params);
  }
  else if (params->plot_type == 1) {
    writeParaview(Solver,params);
  }

}

void writeCSV(solver *Solver, input *params)
{
  ofstream dataFile;
  int iter = params->iter;

  char fileNameC[50];
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
  }

  dataFile.close();
}

void writeParaview(solver *Solver, input *params)
{
  ofstream dataFile;
  int iter = params->iter;

  char fileNameC[50];
  string fileName = params->dataFileName;

#ifndef _NO_MPI
  /* --- All processors write their solution to their own .vtu file --- */
  sprintf(fileNameC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],iter,params->rank);

  /* --- Write 'master' .pvtu file --- */
  if (params->rank == 0) {
    ofstream pVTU;
    char pvtuC[50];
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

    char filnameTmpC[50];
    for (int p=0; p<params->nproc; p++) {
      sprintf(filnameTmpC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],iter,&fileName[0],iter,p);
      pVTU << "    <Piece Source=\"" << string(filnameTmpC) << "\" />" << endl;
    }
    pVTU << "  </PUnstructuredGrid>" << endl;
    pVTU << "</VTKFile>" << endl;

    pVTU.close();

    char datadirC[50];
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
#else
  /* --- Filename to write to --- */
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],iter);
#endif

  dataFile.open(fileNameC);
  dataFile.precision(16);

  if (params->rank == 0)
    cout << "Writing ParaView file " << string(fileNameC) << "...  " << flush;

  // File header
  dataFile << "<?xml version=\"1.0\" ?>" << endl;
  dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
  dataFile << "	<UnstructuredGrid>" << endl;

  // If this is the initial file, need to extrapolate solution to flux points
  if (params->iter==params->initIter) Solver->extrapolateU();

  Solver->extrapolateUMpts();

  if (params->equation == NAVIER_STOKES) {
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
    e.getGridVelPlot(gridVelPpts);
    ppts = e.getPpts();

    if (params->equation == NAVIER_STOKES && params->calcEntropySensor)
      e.getEntropyErrPlot(errPpts);

    int nSubCells = (e.order+2)*(e.order+2);
    int nPpts = (e.order+3)*(e.order+3);
    int nPpts1D = e.order+3;

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
        dataFile << vPpts(k,3) << " ";
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
        dataFile << errPpts(k) << " ";
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

    for (int i=0; i<nPpts1D-1; i++) {
      for (int j=0; j<nPpts1D-1; j++) {
        dataFile << i*nPpts1D     + j   << " ";
        dataFile << i*nPpts1D     + j+1 << " ";
        dataFile << (i+1)*nPpts1D + j+1 << " ";
        dataFile << (i+1)*nPpts1D + j   << " ";
        dataFile << endl;
      }
    }
    dataFile << "				</DataArray>" << endl;

    // Write cell-node offsets
    dataFile << "				<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    for(int k=0; k<nSubCells; k++){
      dataFile << (k+1)*4 << " ";
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    // Write VTK element type
    // 5 = tri, 9 = quad, 10 = tet, 12 = hex
    dataFile << "				<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    for(int k=0; k<nSubCells; k++) {
      dataFile << 9 << " ";
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
    if (iter==1 || (iter/params->monitor_res_freq)%25==0) {
      cout << endl;
      cout << setw(8) << left << "Iter";
      if (params->equation == ADVECTION_DIFFUSION) {
        cout << " Residual " << endl;
      }else if (params->equation == NAVIER_STOKES) {
        cout << setw(colW) << left << "rho";
        cout << setw(colW) << left << "rhoU";
        cout << setw(colW) << left << "rhoV";
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
