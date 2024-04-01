// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef FIELDS_CHAIN_HEADER_H
#define FIELDS_CHAIN_HEADER_H

#include "vector_utils.h"

namespace fields_space{
  // Initializes the couplings for a system of fields on a 1d lattice 
  vec3d get_coupling_matrix_chain(vec1d couplings);

  // Calculates the positions of the nearest neighbours of site r
  vec1i get_neighbours_chain(int r, int Lx);

  // Calculates the direction index of the bond between r1 and r2
  int get_bond_direction_chain(int r1, int r2, int Lx);
}

#endif
