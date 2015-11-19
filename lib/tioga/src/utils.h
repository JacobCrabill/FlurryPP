#pragma once

#include "codetypes.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void kaiser_wrap_(double *,int *,int *,double *,double *,double *,int *);

extern void median_(int *,double *,int *,double *);

/* ---- Oriented Bounding Box Functions ---- */

//!find oriented bounding box for a given set of points
void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);

//! check if a point is inside the provided hole map
int checkHoleMap(double *x,int *nx,int *sam,double *extents);

/*!
 * fill a given hole map using iterative
 * flood fill from outside the marked boundary.
 * boundary is marked by "2"
 */
void fillHoleMap(int *holeMap, int ix[3],int isym);

int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
           double vB[3][3],double xB[3],double dxB[3]);

void writebbox(OBB *obb,int bid);

void writePoints(double *x,int nsearch,int bid);

/* ---- Linked-List Functions ---- */

void deallocateLinkList(DONORLIST *temp);

void deallocateLinkList2(INTEGERLIST *temp);

void insertInList(DONORLIST **donorList,DONORLIST *temp1);

/* ---- Math Functions ---- */

void solvec(double **a,double *b,int *iflag,int n);

void newtonSolve(double f[8][3],double *u1,double *v1,double *w1);

/* ---- Geometry Functions ---- */

void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);

double computeCellVolume(double xv[8][3],int nvert);

#ifdef __cplusplus
}
#endif
