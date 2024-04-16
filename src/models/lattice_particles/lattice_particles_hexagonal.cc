#include "lattice_particles_hexagonal.h"
#include "lattice_particles_geometry.h"
#include "vector_utils.h"

namespace lattice_particles_space {
void get_neighbours_hexagonal(Neighbours &neighbours, int r, int lx, int ly) {
  int i, j;

  array_space::r_to_ij(r, i, j, lx);

  // Left and right coordinates, taking into account lattice edges
  int im = array_space::mod(i - 1, lx);
  int ip = array_space::mod(i + 1, lx);
  int jm = array_space::mod(j - 1, ly);
  int jp = array_space::mod(j + 1, ly);

  array_space::ij_to_r(neighbours[0], ip, j, lx);
  array_space::ij_to_r(neighbours[1], ip, jp, lx);
  array_space::ij_to_r(neighbours[2], i, jp, lx);
  array_space::ij_to_r(neighbours[3], im, j, lx);
  array_space::ij_to_r(neighbours[4], im, jm, lx);
  array_space::ij_to_r(neighbours[5], i, jm, lx);
}

int get_bond_direction_hexagonal(int r1, int r2, int Lx, int Ly) {

  int i1, j1;
  array_space::r_to_ij(r1, i1, j1, Lx);

  int i1m = array_space::mod(i1 - 1, Lx);
  int i1p = array_space::mod(i1 + 1, Lx);
  int j1m = array_space::mod(j1 - 1, Ly);
  int j1p = array_space::mod(j1 + 1, Ly);

  int r1_0, r1_1, r1_2, r1_3, r1_4, r1_5;
  array_space::ij_to_r(r1_0, i1p, j1, Lx);
  array_space::ij_to_r(r1_1, i1p, j1p, Lx);
  array_space::ij_to_r(r1_2, i1, j1p, Lx);
  array_space::ij_to_r(r1_3, i1m, j1, Lx);
  array_space::ij_to_r(r1_4, i1m, j1m, Lx);
  array_space::ij_to_r(r1_5, i1, j1m, Lx);

  if (r2 == r1_0) {
    return 0;
  } else if (r2 == r1_1) {
    return 1;
  } else if (r2 == r1_2) {
    return 2;
  } else if (r2 == r1_3) {
    return 3;
  } else if (r2 == r1_4) {
    return 4;
  } else if (r2 == r1_5) {
    return 5;
  }
  // In case things fail spectacularly
  return -1;
}

int get_conjugate_edge_hexagonal(int edge) {
  if (edge == 3)
    return 6;
  else
    return array_space::mod(edge + 3, 6);
}

} // namespace lattice_particles_space
