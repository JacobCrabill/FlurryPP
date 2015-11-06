/** 
 * Generic Alternating Digital Tree For Search Operations
 */

#pragma once

#include <unordered_set>

// forward declaration for instantiation
class MeshBlock;

#include "MeshBlock.h"

class ADT
{
private :
  
  int ndim;          /** < number of dimensions (usually 3 but can be more) */
  int nelem;         /** < number of elements */
  int *adtIntegers;  /** < integers that define the architecture of the tree */
  double *adtReals;  /** < real numbers that provide the extents of each box */
  double *adtExtents; /** < global extents */
  double *coord;          /** < bounding box of each element */

public :
  ADT() {ndim=6;nelem=0;adtIntegers=NULL;adtReals=NULL;adtExtents=NULL;coord=NULL;}

  ~ADT() 
  {
    if (adtIntegers) free(adtIntegers);
    if (adtReals) free(adtReals);
    if (adtExtents) free(adtExtents);
    adtIntegers=NULL;
    adtReals=NULL;
    adtExtents=NULL;
  }

  void clearData(void)
  {
    if (adtIntegers) free(adtIntegers);
    if (adtReals) free(adtReals);
    if (adtExtents) free(adtExtents);
    adtIntegers=NULL;
    adtReals=NULL;
    adtExtents=NULL;
  }

  void buildADT(int d,int nelements,double *elementBbox);  

  //! Search the ADT for the element containint the point xsearch
  void searchADT_point(MeshBlock *mb,int *cellIndex,double *xsearch);

  //! Search the ADT for all elements overlapping with bounding-box bbox
  void searchADT_box(int *elementList, std::unordered_set<int>& icells, double *bbox);
};

//! Recursively search the input ADT for the MeshBlock cell containing xsearch
void searchIntersections(MeshBlock *mb,int *cellIndex,int *adtIntegers,double *adtReals,double *coord,int level,int node,double *xsearch,int nelem,int ndim);

//! Recursively search the input ADT for all MeshBlock cells intersecting with bbox
void searchBoxIntersections(int *elementList, std::unordered_set<int>& icells, int *adtIntegers, double *adtReals, double *coord, int level, int node, double *bbox, int nelem, int ndim);

void buildADTrecursion(double *coord, double *adtReals, double *adtWork, int *adtIntegers, int *eleIDs, int *adtCount, int side, int parent, int level, int ndim, int nelem, int nav);
