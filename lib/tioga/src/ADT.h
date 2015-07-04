/** 
 * Generic Alternating Digital Tree For Search Operations
 */
// forward declaration for instantiation
class MeshBlock; 
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
  ADT() {ndim=6;nelem=0;adtIntegers=NULL;adtReals=NULL;adtExtents=NULL;coord=NULL;};
  ~ADT() 
    {
      if (adtIntegers) free(adtIntegers);
      if (adtReals) free(adtReals);
      if (adtExtents) free(adtExtents);
      adtIntegers=NULL;
      adtReals=NULL;
      adtExtents=NULL;
    };
  void clearData(void)
    {
      if (adtIntegers) free(adtIntegers);
      if (adtReals) free(adtReals);
      if (adtExtents) free(adtExtents);
      adtIntegers=NULL;
      adtReals=NULL;
      adtExtents=NULL;
    };      
  void buildADT(int d,int nelements,double *elementBbox);  
  void searchADT(MeshBlock *mb,int *cellindx,double *xsearch);
};

