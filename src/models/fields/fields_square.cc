#include "fields_square.h"

namespace fields_space{
  
  vec3d get_coupling_matrix(vec1d couplings){
    
    vec3d coupling_matrix;
    
    // Read the coupling matrix for the bond along the +x axis from the input
    // file. The structure of couplings is assumed to be a flattened version of
    // the coupling_matrix.
    vec2d array_0;
    for(int a=0;a<4;a++){
      vec1d row;
      for(int b=0;b<4;b++){
        row.push_back(couplings[1+a+4*b]);
      }
      array_0.push_back(row);
    }

    coupling_matrix.push_back(array_0);

    for(int k=1;k<4;k++){
      vec2d array_k;
      for(int a=0;a<4;a++){
        vec1d row;
        for(int b=0;b<4;b++){
          int ap = array_space::mod(a-k,4);
          int bp = array_space::mod(a+k,4);
          row.push_back(coupling_matrix[0][ap][bp]);
        }
        array_k.push_back(row);
      }
      coupling_matrix.push_back(array_k);
    }

    return coupling_matrix; 
  }

  vec1i get_neighbours(int r, state_struct &state){
    
    int rm = array_space::mod(r-1,state.Lx);
    int rp = array_space::mod(r+1,state.Lx);

    vec1i neighbours = {rp,rm};

    return neighbours;
  }
}
