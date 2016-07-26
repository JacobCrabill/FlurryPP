/*!
 * \file flurry_interface.cpp
 * \brief Interface functions to Flurry solver
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
#include "flurry_interface.hpp"

Flurry *fr = NULL;

namespace flurry {

#ifndef _NO_MPI
void initialize(MPI_Comm comm_in, char *inputFile, int nGrids, int gridID)
{
  if (!fr) fr = new Flurry(comm_in, nGrids, gridID);

  fr->read_input(inputFile);
}
#else
void initialize(char *input_file)
{
  if (!fr) fr = new Flurry();

  fr->read_input(input_file);
}
#endif

void set_flurry_object(Flurry *_Flurry)
{
  delete fr;

  fr = _Flurry;
}

Flurry* get_flurry_object(void)
{
  return fr;
}

void finalize(void)
{
  //fr->write_wall_time();
  delete fr;
}

/* ---- Data-Acess Functions ---- */

BasicGeo get_basic_geo_data(void)
{
  BasicGeo geo;

  fr->get_basic_geo_data(geo.btag,geo.nnodes,geo.xyz,geo.iblank,geo.nwall,
                           geo.nover,geo.wallNodes,geo.overNodes,geo.nCellTypes,
                           geo.nvert_cell,geo.nCells_type,geo.c2v);

  return geo;
}

double get_q_spt(int ele, int spt, int var)
{
  return fr->get_u_spt(ele,spt,var);
}

double* get_q_spts(void)
{
  return fr->get_u_spts();
}

double* get_q_fpts(void)
{
  return fr->get_u_fpts();
}

ExtraGeo get_extra_geo_data(void)
{
  ExtraGeo geo;

  fr->get_extra_geo_data(geo.nFaceTypes,geo.nvert_face,geo.nFaces_type,
                           geo.f2v,geo.f2c,geo.c2f,geo.iblank_face,
                           geo.iblank_cell,geo.nOverFaces,geo.overFaces,
                           geo.nMpiFaces,geo.mpiFaces,geo.procR,geo.mpiFidR);

  return geo;
}

CallbackFuncs get_callback_funcs(void)
{
  CallbackFuncs call;

  call.get_nodes_per_cell = get_nodes_per_cell;
  call.get_nodes_per_face = get_nodes_per_face;
  call.get_receptor_nodes = get_receptor_nodes;
  call.get_face_nodes = get_face_nodes;
  call.get_q_index_face = get_q_index_face;
  call.donor_inclusion_test = donor_inclusion_test;
  call.donor_frac = donor_frac;
  call.convert_to_modal = convert_to_modal;
  call.get_q_spt = get_q_spt;
  call.get_q_fpt = get_u_fpt;

  return call;
}

/* ---- TIOGA Callback Functions ---- */

void get_nodes_per_cell(int* cellID, int* nNodes)
{
  fr->get_nodes_per_cell(*nNodes);
}

void get_nodes_per_face(int* faceID, int* nNodes)
{
  fr->get_nodes_per_face(*nNodes);
}

void get_receptor_nodes(int* cellID, int* nNodes, double* xyz)
{
  fr->get_receptor_nodes(*cellID, *nNodes, xyz);
}

void get_face_nodes(int* faceID, int* nNodes, double* xyz)
{
  fr->get_face_nodes(*faceID, *nNodes, xyz);
}

void get_q_index_face(int* faceID, int *fpt, int* ind, int* stride)
{
  fr->get_q_index_face(*faceID, *fpt, *ind, *stride);
}

void donor_inclusion_test(int* cellID, double* xyz, int* passFlag, double* rst)
{
  fr->donor_inclusion_test(*cellID, xyz, *passFlag, rst);
}

void donor_frac(int* cellID, double* xyz, int* nweights, int* inode,
                double* weights, double* rst, int* buffsize)
{
  fr->donor_frac(*cellID, *nweights, inode, weights, rst, *buffsize);
}

void convert_to_modal(int* cellID, int* nSpts, double* q_in, int* npts, int* index_out, double* q_out)
{
  //assert(*nSpts == *npts);
  *index_out = (*cellID) * (*nSpts);
  for (int spt = 0; spt < (*nSpts); spt++)
    q_out[spt] = q_in[spt];
}

double &get_u_fpt(int faceID, int fpt, int field)
{
  return fr->get_u_fpt(faceID, fpt, field);
}

} /* namespace flurry */

