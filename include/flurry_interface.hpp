/*!
 * \file flurry_interface.hpp
 * \brief Header file for interface functions to Flurry solver
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014-2016 Jacob Crabill
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
#ifndef _flurry_interface_hpp
#define _flurry_interface_hpp

#include "flurry.hpp"

struct BasicGeo
{
  int btag;         //! Body tag (aka grid ID)
  int nnodes;       //! # of mesh nodes
  double *xyz;      //! Physical positions of all mesh nodes
  int *iblank;      //! Nodal iblank values [to be set externally]
  int nwall;        //! # of wall-boundary nodes
  int nover;        //! # of overset-boundary nodes
  int *wallNodes;   //! List of wall-boundary nodes
  int *overNodes;   //! List of overset-boundary nodes
  int nCellTypes;   //! # of different cell types (hex, tet, prism, etc.) [1 for now]
  int nvert_cell;  //! # of nodes per cell for each cell type
  int nCells_type; //! # of cells for each cell type
  int *c2v;         //! Cell-to-vertex connectivity (one cell type)
};

struct ExtraGeo
{
  int nFaceTypes;   //! # of different face types (quad or tri) [1 for now]
  int nvert_face;  //! # of nodes per face
  int nFaces_type; //! # of faces for each face type
  int *f2v;         //! Face-to-vertex connectivity (one face type)
  int *f2c;         //! Face-to-cell connectivity
  int *c2f;         //! Cell-to-face connectivity
  int *iblank_cell; //! Cell iblank values
  int *iblank_face; //! Face iblank values
  int nOverFaces;   //! # of explicitly-defined overset faces
  int nMpiFaces;    //! # of MPI faces
  int *overFaces;   //! List of explicitly-defined overset faces
  int *mpiFaces;    //! List of MPI face ID's on this rank
  int *procR;       //! Opposite rank for each MPI face
  int *mpiFidR;     //! Face ID of MPI face on opposite rank
};

struct CallbackFuncs
{
  void (*get_nodes_per_cell)(int* cellID, int* nNodes);
  void (*get_nodes_per_face)(int* faceID, int* nNodes);
  void (*get_receptor_nodes)(int* cellID, int* nNodes, double* xyz);
  void (*get_face_nodes)(int* faceID, int* nNodes, double* xyz);
  void (*get_q_index_face)(int* faceID, int *fpt, int* ind, int* stride);
  void (*donor_inclusion_test)(int* cellID, double* xyz, int* passFlag,
                               double* rst);
  void (*donor_frac)(int* cellID, double* xyz, int* nweights, int* inode,
                     double* weights, double* rst, int* buffsize);
  void (*convert_to_modal)(int *cellID, int *nSpts, double *q_in, int *npts,
                           int *index_out, double *q_out);
  double (*get_q_spt)(int cellID, int spt, int var);
  double& (*get_q_fpt)(int faceID, int fpt, int var);
};

namespace flurry {

#ifndef _NO_MPI
void initialize(MPI_Comm comm_in, char* inputFile, int nGrids=1, int gridID=0);
#else
void initialize(char* input_file);
#endif

void set_flurry_object(Flurry *_Flurry);

Flurry* get_flurry_object(void);

void finalize(void);

/* ==== Access functions for mesh data ==== */

BasicGeo get_basic_geo_data(void);

ExtraGeo get_extra_geo_data(void);

CallbackFuncs get_callback_funcs(void);

/* ==== Access functions for solution data ==== */

double get_q_spt(int ele, int spt, int var);
double *get_q_spts(void);
double *get_q_fpts(void);

/* ==== Callback Function Wrappers ==== */

void get_nodes_per_cell(int* cellID, int* nNodes);
void get_nodes_per_face(int* faceID, int* nNodes);
void get_receptor_nodes(int* cellID, int* nNodes, double* xyz);
void get_face_nodes(int* faceID, int* nNodes, double* xyz);
void get_q_index_face(int* faceID, int *fpt, int* ind, int* stride);
void donor_inclusion_test(int* cellID, double* xyz, int* passFlag, double* rst);
void donor_frac(int* cellID, double* xyz, int* nweights, int* inode,
                double* weights, double* rst, int* buffsize);
void convert_to_modal(int *cellID, int *nSpts, double *q_in, int *npts,
                      int *index_out, double *q_out);

double &get_u_fpt(int faceID, int fpt, int field);

} /* namespace Flurry */

#endif // _flurry_interface_hpp

