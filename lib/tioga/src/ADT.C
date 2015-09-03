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

void searchBoxIntersections(MeshBlock *mb,std::set<int> &icells,int *adtIntegers,double *adtReals,
       double *coord,int level,int node,double *bbox,int nelem,int ndim)
{
  int i;
  int d,nodeChild,dimcut;
  double element[ndim];

  for(i=0;i<ndim;i++)
    element[i]=coord[ndim*(adtIntegers[4*node])+i];

  // Check if bbox intersects with bounding box of current mesh element in ADT
  bool flag = true;
  for(i=0;i<3;i++) {
    flag = (flag && (bbox[i+3] >= element[i]  -TOL));
    flag = (flag && (bbox[i]   <= element[i+3]+TOL));
  }

  if (flag) {
    int ind = mb->getCellIndex(adtIntegers[4*node]);
    icells.insert(ind);
  }

  // check the left and right children
  // now
  for(d=1;d<3;d++)
  {
    nodeChild=adtIntegers[4*node+d];
    if (nodeChild > -1) {
      nodeChild = adtIntegers[4*nodeChild+3];

      for(i=0;i<ndim;i++)
        element[i]=adtReals[ndim*nodeChild+i];

      flag = true;
      for(i=0;i<3;i++) {
        flag = (flag && (bbox[i+3] >= element[i]  -TOL));
        flag = (flag && (bbox[i]   <= element[i+3]+TOL));
      }

      if (flag)
        searchBoxIntersections(mb,icells,adtIntegers,adtReals,coord,level+1,
                            nodeChild,bbox,nelem,ndim);
    }
  }
  return;
}

void buildADTrecursion(double *coord,double *adtReals,double *adtWork,int *adtIntegers,
           int *elementsAvailable,int *adtCount,int side,int parent,
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
      adtWork[i]=coord[ndim*elementsAvailable[i]+dimcut];

    // reorder elements with nleft elements to
    // the left of median of adtWork

    median_(elementsAvailable,adtWork,&nav,&coordmid);
    nleft=(nav+1)/2;
    (*adtCount)++;
    ii=(*adtCount)*4;
    adtIntegers[ii]=elementsAvailable[nleft-1];
    adtIntegers[ii+1]=-1;
    adtIntegers[ii+2]=-1;
    adtIntegers[ii+3]=-1;

    // find minimum and maximum bounds of the elements
    // contained in this leaf

    for(i=0;i<nd;i++)
    {
      adtReals[ndim*(*adtCount)+i]=BIGVALUE;
      adtReals[ndim*(*adtCount)+i+nd]=-BIGVALUE;
    }

    for(i=0;i<nav;i++)
    {
      for(j=0;j<nd;j++)
      {
        ii=ndim*(*adtCount)+j;
        iip=ii+nd;
        jj=ndim*elementsAvailable[i]+j;
        jjp=jj+nd;

        adtReals[ii]=min(adtReals[ii],coord[jj]);
        adtReals[iip]=max(adtReals[iip],coord[jjp]);
      }
    }

    // specify that the new element is the child of parent
    // unless root

    if (side > 0)
    {
      adtIntegers[4*parent+side]=elementsAvailable[nleft-1];
    }
    parentToChild=*adtCount;

    // build the left side of the tree

    if (nleft > 1) {
      buildADTrecursion(coord,adtReals,adtWork,adtIntegers,elementsAvailable,
                        adtCount,1,parentToChild,level+1,ndim,nelem,nleft-1);
    }

    // build the right side of the tree

    buildADTrecursion(coord,adtReals,adtWork,adtIntegers,&(elementsAvailable[nleft]),
                      adtCount,2,parentToChild,level+1,ndim,nelem,nav-nleft);
  }
  else if (nav==1) {
    // Base case: only 1 available element left
    (*adtCount)++;
    ii=4*(*adtCount);
    jj=ndim*(*adtCount);
    adtIntegers[ii]=elementsAvailable[0];
    adtIntegers[ii+1]=-1;
    adtIntegers[ii+2]=-1;
    adtIntegers[ii+3]=-1;
    for(j=0;j<ndim;j++)
      adtReals[jj+j]=coord[ndim*elementsAvailable[0]+j];
    if (side > 0) {
      adtIntegers[4*parent+side]=elementsAvailable[0];
    }
  }
}

void ADT::buildADT(int d, int nelements,double *elementBbox)
{
  int i,i2,j6,j,i4;
  int *elementsAvailable;
  double *adtWork;
  int adtCount,parent,level,nav;
  int side;
  double tolerance,delta;

  /* set dimensions and number of elements */

  ndim=d;
  nelem=nelements;

  /* set element bbox pointer */

  coord=elementBbox;

  /*
   * Allocate work arrays
   */

  elementsAvailable=(int *) malloc(sizeof(int)*nelem);
  adtWork=(double *) malloc(sizeof(double)*nelem);

  /*
   * Allocate arrays in the class
   */

  if (adtExtents) free(adtExtents);
  adtExtents=(double *) malloc(sizeof(double)*ndim);
  if (adtIntegers) free(adtIntegers);
  adtIntegers=(int *) malloc(sizeof(int)*4*nelem);
  if (adtReals) free(adtReals);
  adtReals=(double *) malloc(sizeof(double)*nelem*ndim);

  /*
   * Determine extent of elements
   */

  for(i=0;i<ndim/2;i++)
  {
    i2=2*i;
    adtExtents[i2]=BIGVALUE;
    adtExtents[i2+1]=-BIGVALUE;
  }
  for(j=0;j<nelem;j++)
  {
    j6=6*j;
    for(i=0;i<ndim/2;i++)
    {
      i2=2*i;
      adtExtents[i2]=min(adtExtents[i2],coord[j6+i]);
    }
    for(i=0;i<ndim/2;i++)
    {
      i2=2*i+1;
      adtExtents[i2]=max(adtExtents[i2],coord[j6+i+ndim/2]);
    }
  }

  // make the extents 1% larger

  tolerance=0.01;
  for(i=0;i<ndim/2;i++)
  {
    i2=2*i;
    delta=tolerance*(adtExtents[i2+1]-adtExtents[i2]);
    adtExtents[i2]-=delta;
    adtExtents[i2+1]+=delta;
  }

  // Build ADT using a recursive process now

  for(i=0;i<nelem;i++)
    elementsAvailable[i]=i;

  // set initialvalues

  adtCount=-1;
  side=0;
  parent=0;
  level=0;
  nav=nelem;

  buildADTrecursion(coord,adtReals,adtWork,adtIntegers,elementsAvailable,
                    &adtCount,side,parent,level,ndim,nelem,nav);
  //tracei(adtCount);

  // create Inverse map

  //FILE* fp=fopen("adtReals.dat","w");
  //FILE* fp1=fopen("adtInts.dat","w");
  for(i=0;i<nelem;i++)
  {
    i4=4*adtIntegers[4*i];
    adtIntegers[i4+3]=i;
  }

  //for(i=0;i<nelem;i++)
  // {
  //   fprintf(fp,"%.8e %.8e %.8e %.8e %.8e %.8e\n",adtReals[6*i],adtReals[6*i+1],adtReals[6*i+2],adtReals[6*i+3],
   //                                 adtReals[6*i+4],adtReals[6*i+5]);
   // fprintf(fp1,"%d %d %d %d\n",adtIntegers[4*i],adtIntegers[4*i+1],adtIntegers[4*i+2],adtIntegers[4*i+3]);
  // }
  //fclose(fp);
  //fclose(fp1);

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

void ADT::searchADT_box(MeshBlock *mb, std::set<int> &icells, double *bbox)
{
  int rootNode=0;
  icells.clear();

  // Check if the given bounding box intersects with the the bounds of the ADT
  bool flag = true;
  for(int i=0;i<3;i++) {
    flag = (flag && (bbox[i+3] >= adtExtents[2*i]  -TOL));
    flag = (flag && (bbox[i]   <= adtExtents[2*i+1]+TOL));
  }

  // Call recursive routine to check intersections with ADT nodes
  if (flag) searchBoxIntersections(mb,icells,adtIntegers,adtReals,
        coord,0,rootNode,bbox,nelem,ndim);
}

