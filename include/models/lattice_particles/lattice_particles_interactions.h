#ifndef LATTICEPARTICLES_INTERACTIONS_H
#define LATTICEPARTICLES_INTERACTIONS_H

#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"
#include "lattice_particles_geometry.h"

#include <iostream>
#include <cmath>

#include "vector_utils.h"

namespace lattice_particles_space {
// TODO FIgure out, after going thru the code, if this struct is actually
// necessary!
// Structure containing the characteristics of the model interactions:
// couplings          - Interaction energy between neighbouring particles
// energy             - Energy of the system
//
template <int N>
struct interactions_struct {
    ContactMap couplings {};
    double energy {};
    Neighbours<N> neighbours {};
};

// Calculate interactions characteristics of the current state of the system
template <int N>
void initialize_interactions(state_struct &state,
                             interactions_struct<N> &interactions,
                             model_parameters_struct &parameters);

// Print the summary of the interactions characteristics
template <int N>
void print_interactions(state_struct &state, interactions_struct<N> &interactions);

template <int N>
void print_energy(state_struct &state, interactions_struct<N> &interactions);

// Get the contact energy between 2 neighbouring sites, center_site and
// neighbour_site, touching along the edge of center_site contact_edge.
// contact_map is the object containing the energy of all possible contacts
// between pairs of particles.
double get_contact_energy(site_state &center_site, site_state &neighbour_site,
                          int center_site_edge, int neighbour_site_edge,
                          ContactMap contact_map, int n_states);

double get_site_energy(int site_index, state_struct &state,
                       ContactMap coupling_matrix);

template <int N>
double get_energy(state_struct &state, interactions_struct<N> interactions);

} // lattice_particles_space

#endif
