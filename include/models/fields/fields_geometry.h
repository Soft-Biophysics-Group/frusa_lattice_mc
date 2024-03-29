#ifndef FIELDS_GEOMETRY_HEADER_H
#define FIELDS_GEOMETRY_HEADER_H

#include "vector_utils.h"

namespace fields_space{
  // Initializes the couplings for a system of fields on a given lattice
  vec3d get_coupling_matrix(vec1d couplings);

  // Calculates the positions of the nearest neighbours of site r
  vec1i get_neighbours(int r, int Lx, int Ly, int Lz);

  // Calculates the direction index of the bond between r1 and r2
  int get_bond_direction(int r1, int r2, int Lx, int Ly, int Lz);
}

#endif
