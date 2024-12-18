// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "fields_geometry.h"

#if defined FIELDS_CHAIN
#include "fields_chain.h"
#elif defined FIELDS_SQUARE
#include "fields_square.h"
#elif defined FIELDS_HEXAGONAL
#include "fields_hexagonal.h"
#endif

namespace fields_space{
  
  vec3d get_coupling_matrix(vec1d couplings){
#if defined FIELDS_CHAIN
    vec3d coupling_matrix = get_coupling_matrix_chain(couplings);
#elif defined FIELDS_SQUARE
    vec3d coupling_matrix = get_coupling_matrix_square(couplings);
#elif defined FIELDS_HEXAGONAL
    vec3d coupling_matrix = get_coupling_matrix_hexagonal(couplings);
#endif 
    return coupling_matrix;
  }

  vec1i get_neighbours(int r, int Lx, int Ly, int Lz){
#if defined FIELDS_CHAIN
    vec1i neighbours = get_neighbours_chain(r,Lx);
#elif defined FIELDS_SQUARE
    vec1i neighbours = get_neighbours_square(r,Lx,Ly);
#elif defined FIELDS_HEXAGONAL
    vec1i neighbours = get_neighbours_hexagonal(r,Lx,Ly);
#endif 
    return neighbours;
  }
  
  int get_bond_direction(int r1, int r2, int Lx, int Ly, int Lz){
#if defined FIELDS_CHAIN
    int dr = get_bond_direction_chain(r1,r2,Lx);
#elif defined FIELDS_SQUARE
    int dr = get_bond_direction_square(r1,r2,Lx,Ly);
#elif defined FIELDS_HEXAGONAL
    int dr = get_bond_direction_hexagonal(r1,r2,Lx,Ly);
#endif 
    return dr;
  }
}

