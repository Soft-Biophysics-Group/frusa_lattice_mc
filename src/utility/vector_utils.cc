#include "vector_utils.h"

namespace array_space{
  
  void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly, int Lz){
    r = i + Lx*j + Lx*Ly*k;
  }

  void r_to_ijk(int r, int &i, int &j, int &k, int Lx, int Ly, int Lz){
    k = r/(Lx*Ly);
    j = (r-Lx*Ly*k)/Lx;
    i = r-Lx*j-Lx*Ly*k;
  }
}
