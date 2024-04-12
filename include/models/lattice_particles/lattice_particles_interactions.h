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
void initialize_interactions(interactions_struct &interactions,
                             const state_struct &state,
                             const model_parameters_struct &parameters);

// Print the summary of the interactions characteristics
void print_interactions(state_struct &state, interactions_struct &interactions);

void print_energy(state_struct &state, interactions_struct &interactions);

// Get the contact energy between 2 neighbouring sites with indices site1 and
// site2, touching along respective edge indices edge1 and edge2.
// contact_map is the object containing the energy of all possible contacts
// between pairs of particles.
double get_contact_energy(ContactMap contact_map, int state1, int state2,
                          int edge1, int edge2, int n_states);
// overload of previous function for sites on the state
double get_contact_energy(state_struct &state, int site1, int site2,
                          interactions_struct &interactions);

double get_site_energy(state_struct &state, interactions_struct &interactions,
                       int site_index);

double get_energy(state_struct& state, interactions_struct interactions);

} // lattice_particles_space

#endif
