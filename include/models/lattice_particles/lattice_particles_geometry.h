#ifndef LATTICEPARTICLES_GEOMETRY_H
#define LATTICEPARTICLES_GEOMETRY_H

#include "lattice_particles_state.h"
#include "vector_utils.h"

namespace lattice_particles_space {
#if defined LATTICEPARTICLES_HEXAGONAL
using Neighbours = std::array<int, 6>;
#endif

// To be removed once specialized
using Neighbours = std::array<std::size_t, 6>;

// Calculates the positions of the nearest neighbours of site r
void get_neighbours(Neighbours &neighbours, int r, int Lx, int Ly, int Lz);
void get_neighbours(Neighbours &neighbours, std::size_t r, state_struct &state);

// Calculates the direction index of the bond between r1 and r2
int get_bond_direction(int r1, int r2, int Lx, int Ly, int Lz);
int get_bond_direction(std::size_t r1, std::size_t r2, state_struct& state);
} // namespace lattice_particles_space

#endif
