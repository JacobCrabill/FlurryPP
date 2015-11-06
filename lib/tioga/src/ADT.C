/**
 * Build an alternating digital tree
 */
#include <stdio.h>
#include <stdlib.h>
#include "codetypes.h"
#include "ADT.h"

#include "utils.h"

void searchIntersections(MeshBlock *mb,int *cellIndex,int *adtIntegers,double *adtReals,
       double *coord,int level,int node,double *xsearch,int nelem,int ndim)
{
  int i;
  int d,nodeChild,dimcut;
  double element[ndim];
  bool flag;

  for(i=0;i<ndim;i++)
    element[i]=coord[ndim*(adtIntegers[4*node])+i];

  // Check if point is inside bounding box of current mesh element, index 'node'
  flag=1;
  for(i=0;i<ndim/2;i++)
    flag = (flag && (xsearch[i] >=element[i]-TOL));
  for(i=ndim/2;i<ndim;i++)
    flag = (flag && (xsearch[i-ndim/2] <=element[i]+TOL));

  if (flag)
    {
      mb->checkContainment(cellIndex,adtIntegers[4*node],xsearch);
      if (*cellIndex > -1) return;
    }

  // check the left and right children
  // now
  for(d=1;d<3;d++)
  {
    nodeChild=adtIntegers[4*node+d];
    if (nodeChild > -1) {
      nodeChild=adtIntegers[4*nodeChild+3];
      for(i=0;i<ndim;i++)
      {
        element[i]=adtReals[ndim*nodeChild+i];
      }
      flag=1;
      for(i=0;i<ndim/2;i++)
        flag = (flag && (xsearch[i] >=element[i]-TOL));
      for(i=ndim/2;i<ndim;i++)
        flag = (flag && (xsearch[i-ndim/2] <=element[i]+TOL));
      if (flag)
      {
        searchIntersections(mb,cellIndex,adtIntegers,adtReals,coord,level+1,
                            nodeChild,xsearch,nelem,ndim);
        if (*cellIndex > -1) return;
      }
    }
  }
  return;
}

void searchBoxIntersections(int *elementList,std::unordered_set<int> &icells,int *adtIntegers,double *adtReals,
       double *coord,int level,int node,double *bbox,int nelem,int ndim)
{
  double eleBox[ndim];
  for(int i=0;i<ndim;i++)
    eleBox[i] = coord[ndim*(adtIntegers[4*node])+i];

  // Check if bbox intersects with bounding box of current mesh element in ADT
  bool flag = true;
  for (int i=0; i<ndim/2; i++) {
    flag = (flag && (bbox[i+ndim/2] >= eleBox[i]  -TOL));
    flag = (flag && (bbox[i]   <= eleBox[i+ndim/2]+TOL));
  }

  if (flag) {
    int ind = elementList[adtIntegers[4*node]];
    icells.insert(ind);
  }

  // check the left and right children
  // now
  for (int d=1; d<3; d++)
  {
    int nodeChild=adtIntegers[4*node+d];
    if (nodeChild > -1) {
      nodeChild = adtIntegers[4*nodeChild+3];

      double adtBox[ndim];
      for (int i=0; i<ndim; i++)
        adtBox[i]=adtReals[ndim*nodeChild+i];

      flag = true;
      for (int i=0; i<ndim/2; i++) {
        flag = (flag && (bbox[i+ndim/2] >= adtBox[i]  -TOL));
        flag = (flag && (bbox[i] <= adtBox[i+ndim/2]+TOL));
      }

      if (flag) {
        searchBoxIntersections(elementList,icells,adtIntegers,adtReals,coord,level+1,
                               nodeChild,bbox,nelem,ndim);
      }
    }
  }
}

void buildADTrecursion(double *coord,double *adtReals,double *adtWork,int *adtIntegers,
           int *eleIDs,int *adtCount,int side,int parent,
           int level,int ndim,int nelem, int nav)
{

  /* RECALL:
   * --coord contains the bounding box for each individual *mesh* element
   * --adtReals contains the bounding box for each *ADT* element/node
   */
  int nd=ndim/2;
  double coordmid;
  int i,j;
  int dimcut;
  int nleft;
  int ii,iip,jj,jjp;
  int parentToChild;

  if (nav > 1) {

    // find the dimension to create the cut

    dimcut=(level%ndim);

    // collect coordinates along the dimension dimcut

    for(i=0;i<nav;i++)
      adtWork[i]=coord[ndim*eleIDs[i]+dimcut];

    // reorder elements with nleft elements to
    // the left of median of adtWork

    median_(eleIDs,adtWork,&nav,&coordmid);
    nleft=(nav+1)/2;
    (*adtCount)++;
    ii=(*adtCount)*4;
    adtIntegers[ii]=eleIDs[nleft-1];
    adtIntegers[ii+1]=-1;
    adtIntegers[ii+2]=-1;
    adtIntegers[ii+3]=-1;

    // find minimum and maximum bounds of the elements
    // contained in this leaf

    for(i=0;i<nd;i++)
    {
      adtReals[ndim*(*adtCount)+i]   = BIGVALUE;
      adtReals[ndim*(*adtCount)+i+nd]=-BIGVALUE;
    }

    for(i=0;i<nav;i++)
    {
      for(j=0;j<nd;j++)
      {
        ii=ndim*(*adtCount)+j;
        iip=ii+nd;
        jj=ndim*eleIDs[i]+j;
        jjp=jj+nd;

        adtReals[ii]=min(adtReals[ii],coord[jj]);
        adtReals[iip]=max(adtReals[iip],coord[jjp]);
      }
    }

    // specify that the new element is the child of parent
    // unless root

    if (side > 0)
    {
      adtIntegers[4*parent+side]=eleIDs[nleft-1];
    }
    parentToChild=*adtCount;

    // build the left side of the tree

    if (nleft > 1) {
      buildADTrecursion(coord,adtReals,adtWork,adtIntegers,eleIDs,
                        adtCount,1,parentToChild,level+1,ndim,nelem,nleft-1);
    }

    // build the right side of the tree

    buildADTrecursion(coord,adtReals,adtWork,adtIntegers,&(eleIDs[nleft]),
                      adtCount,2,parentToChild,level+1,ndim,nelem,nav-nleft);
  }
  else if (nav==1) {
    // Base case: only 1 available element left
    (*adtCount)++;
    ii=4*(*adtCount);
    jj=ndim*(*adtCount);
    adtIntegers[ii]=eleIDs[0];
    adtIntegers[ii+1]=-1;
    adtIntegers[ii+2]=-1;
    adtIntegers[ii+3]=-1;
    for(j=0;j<ndim;j++)
      adtReals[jj+j]=coord[ndim*eleIDs[0]+j];
    if (side > 0) {
      adtIntegers[4*parent+side]=eleIDs[0];
    }
  }
}

void ADT::buildADT(int d, int nelements,double *elementBbox)
{
  /* set dimensions and number of elements */

  ndim=d;
  nelem=nelements;

  /* set element bbox pointer */

  coord=elementBbox;

  /* Allocate work arrays */

  int *elementsAvailable=(int *) malloc(sizeof(int)*nelem);
  double *adtWork=(double *) malloc(sizeof(double)*nelem);

  /* Allocate arrays in the class */

  if (adtExtents) free(adtExtents);
  adtExtents=(double *) malloc(sizeof(double)*ndim);
  if (adtIntegers) free(adtIntegers);
  adtIntegers=(int *) malloc(sizeof(int)*4*nelem);
  if (adtReals) free(adtReals);
  adtReals=(double *) malloc(sizeof(double)*nelem*ndim);

  /* Determine extent of elements */

  for(int i=0; i<ndim/2; i++)
  {
    int i2=2*i;
    adtExtents[i2]=BIGVALUE;
    adtExtents[i2+1]=-BIGVALUE;
  }

  for(int j=0; j<nelem; j++)
  {
    int jd=ndim*j;
    for(int i=0; i<ndim/2; i++)
    {
      int i2=2*i;
      adtExtents[i2]=min(adtExtents[i2],coord[jd+i]);
    }
    for(int i=0; i<ndim/2; i++)
    {
      int i2=2*i+1;
      adtExtents[i2]=max(adtExtents[i2],coord[jd+i+ndim/2]);
    }
  }

  // make the extents 1% larger

  double tolerance=0.01;
  for(int i=0; i<ndim/2; i++)
  {
    int i2=2*i;
    double delta=tolerance*(adtExtents[i2+1]-adtExtents[i2]);
    adtExtents[i2]-=delta;
    adtExtents[i2+1]+=delta;
  }

  // Build ADT using a recursive process now

  for(int i=0; i<nelem; i++)
    elementsAvailable[i]=i;

  // set initialvalues

  int adtCount=-1;
  int side=0;
  int parent=0;
  int level=0;
  int nav=nelem;

  buildADTrecursion(coord,adtReals,adtWork,adtIntegers,elementsAvailable,
                    &adtCount,side,parent,level,ndim,nelem,nav);

  // create Inverse map [ADT index <== original ele ID]
  // adtInt[eleID+3] = adtInd
  for(int i=0; i<nelem; i++)
  {
    int eleID = 4*adtIntegers[4*i];
    adtIntegers[eleID+3] = i;
  }

  free(elementsAvailable);
  free(adtWork);
}


void ADT::searchADT_point(MeshBlock *mb, int* cellIndex, double *xsearch)
{
  int rootNode=0;
  *cellIndex=-1;

  // check if the given point is in the bounds of the ADT
  int flag=1;
  for(int i=0;i<ndim/2;i++)
    flag = (flag && (xsearch[i] >= adtExtents[2*i]-TOL));
  for(int i=0;i<ndim/2;i++)
    flag= (flag && (xsearch[i] <= adtExtents[2*i+1]+TOL));

  // call recursive routine to check intersections with ADT nodes
  if (flag) searchIntersections(mb,cellIndex,adtIntegers,adtReals,
        coord,0,rootNode,xsearch,nelem,ndim);
}

void ADT::searchADT_box(int *elementList, std::unordered_set<int> &icells, double *bbox)
{
  int rootNode=0;
  icells.clear();

  // Check if the given bounding box intersects with the the bounds of the ADT
  bool flag = true;
  for(int i=0;i<ndim/2;i++) {
    flag = (flag && (bbox[i+ndim/2] >= adtExtents[2*i]  -TOL));
    flag = (flag && (bbox[i]   <= adtExtents[2*i+1]+TOL));
  }

  // Call recursive routine to check intersections with ADT nodes
  if (flag) searchBoxIntersections(elementList,icells,adtIntegers,adtReals,
        coord,0,rootNode,bbox,nelem,ndim);
}

void ADT::searchADT_box_debug(int *elementList, std::unordered_set<int> &icells, double *bbox)
{
  int rootNode=0;
  icells.clear();

  cout << "ADT: Search box: " << bbox[0] << ", " << bbox[1] << ", " << bbox[2] << ", " << bbox[3] << " | ndim = " << ndim << endl;

  for (int i=0; i<nelem; i++) {
    cout << i << ": " << coord[ndim*i+0] << ", " << coord[ndim*i+1] << ", " << coord[ndim*i+2] << ", " << coord[ndim*i+3] << endl;
  }
  for (int i=0; i<nelem; i++) {
    cout << "adt i " << i << ": ele " << adtIntegers[4*i] << endl;
  }

  // Check if the given bounding box intersects with the the bounds of the ADT
  bool flag = true;
  for(int i=0;i<ndim/2;i++) {
    flag = (flag && (bbox[i+ndim/2] >= adtExtents[2*i]  -TOL));
    flag = (flag && (bbox[i]   <= adtExtents[2*i+1]+TOL));
  }

  // Call recursive routine to check intersections with ADT nodes
  if (flag) searchBoxIntersections_debug(elementList,icells,adtIntegers,adtReals,
        coord,0,rootNode,bbox,nelem,ndim);
}

void searchBoxIntersections_debug(int *elementList,std::unordered_set<int> &icells,int *adtIntegers,double *adtReals,
       double *coord,int level,int node,double *bbox,int nelem,int ndim)
{
  double eleBox[ndim];
  for(int i=0;i<ndim;i++)
    eleBox[i] = coord[ndim*(adtIntegers[4*node])+i];

  cout << "ADT bbox for ele " << adtIntegers[4*node] << ": ";
  for (int i=0; i<ndim; i++)
    cout << eleBox[i] << ", ";
  cout << endl;
  int e1 = adtIntegers[4*node+1];
  int e2 = adtIntegers[4*node+2];
  int a1 = adtIntegers[4*e1+3];
  int a2 = adtIntegers[4*e2+3];
  cout << "  Children: ele IDs " << e1 << ", " << e2 << "; adt IDs " << a1 << ", " << a2 << endl;

  // Check if bbox intersects with bounding box of current mesh element in ADT
  bool flag = true;
  for (int i=0; i<ndim/2; i++) {
    flag = (flag && (bbox[i+ndim/2] >= eleBox[i]  -TOL));
    flag = (flag && (bbox[i]   <= eleBox[i+ndim/2]+TOL));
  }

  if (flag) {
    int ind = elementList[adtIntegers[4*node]];
    icells.insert(ind);
    cout << "  Found in ADT: ind " << ind << endl;
  }

  // check the left and right children
  // now
  for (int d=1; d<3; d++)
  {
    int nodeChild=adtIntegers[4*node+d];
    if (nodeChild > -1) {
      nodeChild = adtIntegers[4*nodeChild+3];

      double adtBox[ndim];
      for (int i=0; i<ndim; i++)
        adtBox[i]=adtReals[ndim*nodeChild+i];

      flag = true;
      for (int i=0; i<ndim/2; i++) {
        flag = (flag && (bbox[i+ndim/2] >= adtBox[i]  -TOL));
        flag = (flag && (bbox[i] <= adtBox[i+ndim/2]+TOL));
      }

      if (flag) {
        searchBoxIntersections_debug(elementList,icells,adtIntegers,adtReals,coord,level+1,
                               nodeChild,bbox,nelem,ndim);
      }
    }
  }
}
