
#include "global.hpp"

class mesh
{
public:
  //! Default constructor
  mesh();

  void setup();

  //! Default destructor
  ~mesh();

  void get_loc_spts(int eType, int order);

  void get_loc_fpts(int eType, int order);

private:
  string spts_type_tri, fpts_type_tri;
  string spts_type_quad, fpts_type_quad;

};
