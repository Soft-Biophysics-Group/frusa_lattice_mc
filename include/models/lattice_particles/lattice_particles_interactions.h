#ifndef LATTICEPARTICLES_INTERACTIONS_H
#define LATTICEPARTICLES_INTERACTIONS_H

#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"
#include "lattice_particles_geometry.h"

#include <iostream>
#include <cmath>
#include <array>

#include "vector_utils.h"

namespace lattice_particles_space {
// TODO FIgure out, after going thru the code, if this struct is actually
// necessary!
// Structure containing the characteristics of the model interactions:
// couplings          - Interaction energy between neighbouring particles
// energy             - Energy of the system

struct interactions_struct {
    ContactMap couplings {};
    double energy {};
    Neighbours neighbours {};
    int n_neighbours {};
};

// Calculate interactions characteristics of the current state of the system
void initialize_interactions(state_struct &state,
                             interactions_struct &interactions,
                             model_parameters_struct &parameters);

// Print the summary of the interactions characteristics
void print_interactions(state_struct &state, interactions_struct &interactions);

void print_energy(state_struct &state, interactions_struct &interactions);

// Get the contact energy between 2 neighbouring sites, center_site and
// neighbour_site, touching along the edge of center_site contact_edge.
// contact_map is the object containing the energy of all possible contacts
// between pairs of particles.
double get_contact_energy(site_state &center_site, site_state &neighbour_site,
                          int center_site_edge, int neighbour_site_edge,
                          ContactMap contact_map, int n_states);
double get_contact_energy(std::size_t center_index, std::size_t neighbour_index,
                          state_struct& state, interactions_struct& interactions);

double get_site_energy(site_state &site, state_struct &state,
                       interactions_struct& interactions);

double get_energy(state_struct& state, interactions_struct interactions);

} // lattice_particles_space

#endif
