#include <array>
#include "lattice_particles_geometry.h"

#if defined LATTICEPARTICLES_HEXAGONAL
#include "lattice_particles_hexagonal.h"
#endif

namespace lattice_particles_space {

template <int N>
void get_neighbours(Neighbours<N>& neighbours, int r, int Lx, int Ly, int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
    get_neighbours_hexagonal(neighbours, r, Lx, Ly);
#endif
}

int get_bond_directions(int r1, int r2, int Lx, int Ly, int Lz) {
#if defined LATTICEPARTICLES_HEXAGONAL
        return get_bond_directions_hexagonal(int r1, int r2, int Lx, int Ly);
#endif
    return dr;
}
} // namespace lattice_particles_space
