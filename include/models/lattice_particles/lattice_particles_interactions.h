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
    int n_edges {};
};

// Calculate interactions characteristics of the current state of the system
void initialize_interactions(interactions_struct &interactions,
                             state_struct &state,
                             model_parameters_struct &parameters);

// Print the summary of the interactions characteristics
void print_interactions(state_struct &state, interactions_struct &interactions);

void print_energy(state_struct &state, interactions_struct &interactions);

// Get the contact energy between 2 neighbouring sites with indices site1 and
// site2, touching along respective edge indices edge1 and edge2.
// contact_map is the object containing the energy of all possible contacts
// between pairs of particles.
double get_contact_energy(state_struct &state, int site1, int site2,
                          interactions_struct &interactions);

double get_site_energy(state_struct &state, interactions_struct &interactions,
                       int site_index);

double get_energy(state_struct& state, interactions_struct& interactions);

// This function takes as the reference edge for the contact the smallest
// edge and assigns a contact identity based on it.
int get_contact_index(state_struct &state, int site1, int site2, int edge1,
                      interactions_struct& interactions);

int get_contact_index_full_empty(int edge, int type, int orientation,
                                 int n_edges);
int get_contact_index_full_full(int state1, int state2, int edge1, int n_states,
                                int n_edges, int n_types);

std::ostream& operator<< (std::ostream& out, interactions_struct& interactions); 
} // lattice_particles_space

#endif
