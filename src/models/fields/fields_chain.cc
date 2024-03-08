#include "fields_chain.h"

namespace fields_space{
  
  vec3d get_coupling_matrix(vec1d couplings){
    
    double k11 = couplings[1];
    double k12 = couplings[2];
    double k21 = couplings[3];
    
    vec3d coupling_matrix = {{{k11,k12},{k21,k11}},
                             {{k11,k21},{k12,k11}}};

    return coupling_matrix; 
  }

  vec1i get_neighbours(int r, state_struct &state){
    
    int rm = (state.Lx+((r-1)%state.Lx))%state.Lx;
    int rp = (r+1)%state.Lx;

    vec1i neighbours = {rp,rm};

    return neighbours;
  }
}
