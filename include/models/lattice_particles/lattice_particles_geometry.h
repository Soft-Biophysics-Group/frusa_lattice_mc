#ifndef LATTICEPARTICLES_GEOMETRY_H
#define LATTICEPARTICLES_GEOMETRY_H

#include "vector_utils.h"


namespace lattice_particles_space{
  template <int N>
    using Neighbours = std::array<int, N>;

  // Calculates the positions of the nearest neighbours of site r
  vec1i get_neighbours(int r, int Lx, int Ly, int Lz);

  // Calculates the direction index of the bond between r1 and r2
  int get_bond_direction(int r1, int r2, int Lx, int Ly, int Lz);
} //namespace lattice_particles_space

#endif
