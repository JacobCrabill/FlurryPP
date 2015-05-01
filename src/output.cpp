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
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],iter);

  dataFile.precision(15);
  dataFile.setf(ios_base::fixed);

  dataFile.open(fileNameC);

  // File header
  dataFile << "<?xml version=\"1.0\" ?>" << endl;
  //dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\"" << endl;
  dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
  dataFile << "	<UnstructuredGrid>" << endl;

  for (auto& e:Solver->eles) {
    // If this is the initial file, need to extrapolate solution to flux points
    if (params->iter==0) Solver->extrapolateU();

    Solver->extrapolateUMpts();

    // The combination of spts + fpts will be the plot points
    matrix<double> vPpts;
    vector<point> ppts;
    e.getPrimitivesPlot(vPpts);
    e.getPpts(ppts);

    int nSubCells = (e.order+2)*(e.order+2);
    int nPpts = (e.order+3)*(e.order+3);
    int nPpts1D = e.order+3;

    // Write cell header
    dataFile << "		<Piece NumberOfPoints=\"" << nPpts << "\" NumberOfCells=\"" << nSubCells << "\">" << endl;

    /* ==== Write out solution to file ==== */
    dataFile << "			<PointData>" << endl;

    /* --- Density --- */
    dataFile << "				<DataArray type= \"Float32\" Name=\"Density\" format=\"ascii\">" << endl;
    for(int k=0; k<nPpts; k++) {
      dataFile << vPpts(k,0) << " ";
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    /* --- Pressure --- */
    dataFile << "				<DataArray type= \"Float32\" Name=\"Pressure\" format=\"ascii\">" << endl;
    for(int k=0; k<nPpts; k++) {
      dataFile << vPpts(k,3) << " ";
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    /* --- Velocity --- */
    dataFile << "				<DataArray type= \"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">" << endl;
    for(int k=0; k<nPpts; k++) {
      // Divide momentum components by density to obtain velocity components
      dataFile << vPpts(k,1) << " " << vPpts(k,2) << " ";

      // In 2D the z-component of velocity is not stored, but Paraview needs it so write a 0.
      if(params->nDims==2) {
        dataFile << 0.0 << " ";
      }
      else {
        dataFile << vPpts[k][3] << " ";
      }
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

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

    // Write cell IDs...?
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
}

void writeParaviewBinary(solver *Solver, input *params)
{
  ofstream dataFile;
  int iter = params->iter;

  char fileNameC[50];
  string fileName = params->dataFileName;
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],iter);

  //dataFile.precision(15);
  //dataFile.setf(ios_base::fixed);
  //dataFile.open(fileNameC, std::ofstream::app | std::ofstream::binary);
  dataFile.open(fileNameC, std::ofstream::binary);

  char bufi[sizeof(int)];
  char bufu[sizeof(uint)];
  char bufd[sizeof(double)];

  // File header
  dataFile << "<?xml version=\"1.0\" ?>" << endl;
  dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl;
  //dataFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
  dataFile << "	<UnstructuredGrid>" << endl;

  for (auto& e:Solver->eles) {
    // If this is the initial file, need to extrapolate solution to flux points
    if (params->iter==0) Solver->extrapolateU();

    Solver->extrapolateUMpts();

    // The combination of spts + fpts will be the plot points
    matrix<double> vPpts;
    vector<point> ppts;
    e.getPrimitivesPlot(vPpts);
    e.getPpts(ppts);

    int nSubCells = (e.order+2)*(e.order+2);
    int nPpts = (e.order+3)*(e.order+3);
    int nPpts1D = e.order+3;

    // Write cell header
    dataFile << "		<Piece NumberOfPoints=\"" << nPpts << "\" NumberOfCells=\"" << nSubCells << "\">" << endl;

    /* ==== Write out solution to file ==== */
    dataFile << "			<PointData>" << endl;

    /* --- Density --- */
    dataFile << "				<DataArray type= \"Float32\" Name=\"Density\" format=\"binary\">" << endl;
    for(int k=0; k<nPpts; k++) {
      bitswap(vPpts(k,0),bufd);
      dataFile.write(bufd,sizeof(double));
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    /* --- Pressure --- */
    dataFile << "				<DataArray type= \"Float32\" Name=\"Pressure\" format=\"binary\">" << endl;
    for(int k=0; k<nPpts; k++) {
      bitswap(vPpts(k,3),bufd);
      dataFile.write(bufd,sizeof(double));
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    /* --- Velocity --- */
    dataFile << "				<DataArray type= \"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"binary\">" << endl;
    for(int k=0; k<nPpts; k++) {
      // Divide momentum components by density to obtain velocity components
      bitswap(vPpts(k,1),bufd);
      dataFile.write(bufd,sizeof(double));
      bitswap(vPpts(k,2),bufd);
      dataFile.write(bufd,sizeof(double));

      // In 2D the z-component of velocity is not stored, but Paraview needs it so write a 0.
      if(params->nDims==2) {
        bitswap((double)0.,bufd);
        dataFile.write(bufd,sizeof(double));
      }
      else {
        bitswap(vPpts(k,3),bufd);
        dataFile.write(bufd,sizeof(double));
      }
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    dataFile << "			</PointData>" << endl;

    /* ==== Write Out Cell Points & Connectivity==== */

    /* --- Write out the plot point coordinates --- */
    dataFile << "			<Points>" << endl;
    dataFile << "				<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << endl;

    // Loop over plot points in element
    for(int k=0; k<nPpts; k++) {
      for(int l=0;l<params->nDims;l++) {
        bitswap(ppts[k][l],bufd);
        dataFile.write(bufd,sizeof(double));
      }

      // If 2D, write a 0 as the z-component
      if(params->nDims == 2) {
        bitswap((double)0.,bufd);
        dataFile.write(bufd,sizeof(double));
      }
    }

    dataFile << endl;
    dataFile << "				</DataArray>" << endl;
    dataFile << "			</Points>" << endl;

    /* --- Write out Cell data: connectivity, offsets, element types --- */
    dataFile << "			<Cells>" << endl;

    /* --- Write connectivity array --- */
    dataFile << "				<DataArray type=\"Int16\" Name=\"connectivity\" format=\"binary\">" << endl;

    for (int i=0; i<nPpts1D-1; i++) {
      for (int j=0; j<nPpts1D-1; j++) {
        bitswap((int)(i*nPpts1D + j),bufi);
        dataFile.write(bufi,sizeof(int));
        bitswap((int)(i*nPpts1D + j+1),bufi);
        dataFile.write(bufi,sizeof(int));
        bitswap((int)((i+1)*nPpts1D + j+1),bufi);
        dataFile.write(bufi,sizeof(int));
        bitswap((int)((i+1)*nPpts1D + j),bufi);
        dataFile.write(bufi,sizeof(int));
      }
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    // Write cell IDs...?
    dataFile << "				<DataArray type=\"Int16\" Name=\"offsets\" format=\"binary\">" << endl;
    for(int k=0; k<nSubCells; k++){
      bitswap((int)(k+1)*4,bufi);
      dataFile.write(bufi,sizeof(int));
    }
    dataFile << endl;
    dataFile << "				</DataArray>" << endl;

    // Write VTK element type
    // 5 = tri, 9 = quad, 10 = tet, 12 = hex
    dataFile << "				<DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << endl;
    for(int k=0; k<nSubCells; k++) {
      bitswap((uint)9,bufu);
      dataFile.write(bufu,sizeof(uint));
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
}


void writeResidual(solver *Solver, input *params)
{
  vector<double> res(params->nFields), resTmp(params->nFields);
  int iter = params->iter;

  for (auto& e:Solver->eles) {
    resTmp = e.getResidual(params->resType);
    if(checkNaN(resTmp)) FatalError("NaN Encountered in Solution Residual!");

    for (int i=0; i<params->nFields; i++) {
      if (params->resType == 3) {
        // Infinity norm [max residual over all spts]
        res[i] = max(res[i],resTmp[i]);
      }else{
        res[i] += resTmp[i];
      }
    }
  }

  // If taking 2-norm, res is sum squared; take sqrt to complete
  if (params->resType == 2) {
    for (auto& i:res) i = sqrt(i);
  }

  int colW = 16;
  cout.precision(8);
  cout.setf(ios::fixed, ios::floatfield);
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
    }
    cout << endl;
  }

  cout << setw(8) << left << iter;
  for (int i=0; i<params->nFields; i++) {
    cout << setw(colW) << left << res[i];
  }
  cout << endl;
}
