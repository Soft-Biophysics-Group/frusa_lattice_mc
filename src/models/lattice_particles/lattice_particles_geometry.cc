#include "lattice_particles_geometry.h"
#include <array>

#if defined LATTICEPARTICLES_HEXAGONAL
#include "lattice_particles_hexagonal.h"
#endif

// To be removed once testing is done
#include "lattice_particles_hexagonal.h"

namespace lattice_particles_space {

void get_neighbours(Neighbours &neighbours, int r, int Lx, int Ly,
                    [[maybe_unused]] int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
  get_neighbours_hexagonal(neighbours, r, Lx, Ly);
#endif
  get_neighbours_hexagonal(neighbours, r, Lx, Ly);
}
void get_neighbours(Neighbours &neighbours, int r, state_struct &state) {
  get_neighbours(neighbours, r, state.lx, state.ly, state.lz);
}

int get_bond_direction(int r1, int r2, int Lx, int Ly,
                       [[maybe_unused]] int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
  return get_bond_direction_hexagonal(r1, r2, Lx, Ly);
#endif
  return get_bond_direction_hexagonal(r1, r2, Lx, Ly);
}
int get_bond_direction(int r1, int r2, state_struct &state) {
  return get_bond_direction(static_cast<int>(r1), static_cast<int>(r2),
                            state.lx, state.ly, state.lz);
}

int get_conjugate_edge(int edge) {
#if defined LATTICEPARTICLES_HEXAGONAL
  return get_conjugate_edge_hexagonal(edge);
#endif
  return get_conjugate_edge_hexagonal(edge);
}
} // namespace lattice_particles_space
