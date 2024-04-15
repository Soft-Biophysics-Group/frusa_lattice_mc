// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#include "fields_chain.h"

namespace fields_space{
  
  vec3d get_coupling_matrix_chain(vec1d couplings){
    
    double k11 = couplings[1];
    double k12 = couplings[2];
    double k21 = couplings[3];
    
    vec3d coupling_matrix = {{{k11,k12},{k21,k11}},
                             {{k11,k21},{k12,k11}}};

    return coupling_matrix; 
  }

  vec1i get_neighbours_chain(int r, int Lx){
    
    int rm = array_space::mod(r-1,Lx);
    int rp = array_space::mod(r+1,Lx);

    vec1i neighbours = {rp,rm};

    return neighbours;
  }

  int get_bond_direction_chain(int r1, int r2, int Lx){

    int dr;

    int r1m = array_space::mod(r1-1,Lx);
    int r1p = array_space::mod(r1+1,Lx);
    
    if(r2==r1p){
      dr = 0;
    }
    else{
      dr = 1;
    }
    return dr;
  }
}
