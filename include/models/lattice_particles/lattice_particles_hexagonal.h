#ifndef LATTICEPARTICLES_HEXAGONAL_H
#define LATTICEPARTICLES_HEXAGONAL_H

#include "lattice_particles_geometry.h"
#include "vector_utils.h"

namespace lattice_particles_space {
  void get_neighbours_hexagonal(Neighbours<6>& neighbours, int r, int lx, int ly);
  int get_bond_direction_hexagonal(int r1, int r2, int lx, int ly);
} // namespace lattice_particles_space

#endif
