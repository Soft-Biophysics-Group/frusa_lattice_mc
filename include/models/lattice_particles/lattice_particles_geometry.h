#ifndef LATTICEPARTICLES_GEOMETRY_H
#define LATTICEPARTICLES_GEOMETRY_H

#include "lattice_particles_state.h"
#include "vector_utils.h"

namespace lattice_particles_space {

// To be removed once specialized
// TODO Figure out how to have a constexpr global constant number of edges
// with conditional inclusion
inline constexpr int n_edges{6};
using Neighbours = std::array<int, n_edges>;

// Calculates the positions of the nearest neighbours of site r by updating the
// array Neighbours
void get_neighbours(Neighbours &neighbours, int r, int Lx, int Ly, int Lz);
void get_neighbours(Neighbours &neighbours, int r, state_struct &state);

// Returns the direction index of the bond between r1 and r2
int get_bond_direction(int r1, int r2, int Lx, int Ly, int Lz);
int get_bond_direction(int r1, int r2, state_struct& state);

// Gets the edge facing a given edge
int get_conjugate_edge(int edge);
} // namespace lattice_particles_space

#endif
