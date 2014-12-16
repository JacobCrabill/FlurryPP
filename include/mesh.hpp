
#include "global.hpp"

class mesh
{
public:
  //! Default constructor
  mesh();

  mesh.setup();

  //! Default destructor
  ~mesh();

  get_loc_spts(int eType, int order);

  get_loc_fpts(int eType, int order);

private:
  string spts_type_tri, fpts_type_tri;
  string spts_type_quad, fpts_type_quad;

};
