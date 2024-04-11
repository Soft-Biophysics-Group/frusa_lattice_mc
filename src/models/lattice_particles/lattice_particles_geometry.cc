#include <array>
#include "lattice_particles_geometry.h"

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
  get_neighbours(neighbours, r, state);
}

int get_bond_directions(int r1, int r2, int Lx, int Ly,
                        [[maybe_unused]] int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
        return get_bond_directions_hexagonal(r1, r2, Lx, Ly);
#endif
return get_bond_direction_hexagonal(r1, r2, Lx, Ly);
}
} // namespace lattice_particles_space
