#include <array>
#include "lattice_particles_geometry.h"

#if defined LATTICEPARTICLES_HEXAGONAL
#include "lattice_particles_hexagonal.h"
#endif

namespace lattice_particles_space {

void get_neighbours(Neighbours& neighbours, int r, int Lx, int Ly, int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
    get_neighbours_hexagonal(neighbours, r, Lx, Ly);
#endif
}

int get_bond_directions(int r1, int r2, int Lx, int Ly, int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
        return get_bond_directions_hexagonal(int r1, int r2, int Lx, int Ly);
#endif
}
} // namespace lattice_particles_space
