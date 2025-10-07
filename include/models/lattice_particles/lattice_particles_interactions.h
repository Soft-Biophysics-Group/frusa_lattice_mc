#ifndef LATTICEPARTICLES_INTERACTIONS_H
#define LATTICEPARTICLES_INTERACTIONS_H

#include "lattice_particles_geometry.h"
#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"

#include <array>
#include <cmath>
#include <iostream>

#include "vector_utils.h"

namespace lattice_particles_space {

// Structure containing the characteristics of the model interactions:
// couplings          - Interaction energy between neighbouring particles
// energy             - Total energy of the system
// neighbours         - Array of indices neighbouring a given site, updated at
//                      each energy calculation. Useful to avoid constant
//                      allocations.
// n_edges:           - number of neighbours of lattice sites; hopefully
//                      something I can remove in the future.
struct interactions_struct {
  ContactMap couplings{};
  double energy{};
  Neighbours neighbours{};
  int n_edges{};
};

void initialize_interactions(state_struct &state,
                             interactions_struct &interactions,
                             model_parameters_struct &parameters);

// Print the summary of the interactions characteristics
void print_interactions(state_struct &state, interactions_struct &interactions);

// Print current system energy
void print_energy(state_struct &state, interactions_struct &interactions);

// Get the contact energy between 2 neighbouring sites with indices site1 and
// site2, touching along respective edge indices edge1 and edge2.
double get_contact_energy(state_struct &state, int site1, int site2,
                          interactions_struct &interactions);

// Get total energy of a given site, which is the sum of contact energy with
// its neighbours
double get_site_energy(state_struct &state, interactions_struct &interactions,
                       int site_index);

// Get total energy of the system
double get_energy(state_struct &state, interactions_struct &interactions);

// Returns the contact index for two sites, which is the coefficient of the
// corresponding contact in the coupling vector
int get_contact_index(state_struct &state, int site1, int site2,
                      interactions_struct &interactions);

int get_contact_index_full_empty(int edge, int type, int orientation,
                                 int n_edges);

int get_contact_index_full_full(int state1, int state2, int edge1, int n_states,
                                int n_edges, int n_types);

std::ostream &operator<<(std::ostream &out, interactions_struct &interactions);
} // namespace lattice_particles_space

#endif
