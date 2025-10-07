#ifndef PARTICLES_INTERACTIONS_H
#define PARTICLES_INTERACTIONS_H

/**
* Vincent Ouazan-Reboul, 2024
* Classes and functions to measure, and keep track of, the
* orientation-dependent interactions between particles living on a Bravais
* lattice.
*/

#include "particles_parameters.h"
#include "particles_state.h"

#include <array>
#include <cmath>
#include <iostream>

#include "vector_utils.h"
#include "geometry.h"

namespace particles_space {

// Structure containing the characteristics of the model interactions:
struct interactions_struct {
  // Interaction energy between neighbouring particles, flattened into 1
  // dimension
  vec1d couplings{};
  // Current total energy of the system
  double energy{};
  // Number of neighbours of each lattice site. Imposed by choice of lattice
  int n_edges{};
};

void initialize_interactions(state_struct& state,
                             interactions_struct& interactions,
                             model_parameters_struct& parameters,
                             geometry_space::Geometry& geometry);

// Printers
void print_interactions(interactions_struct &interactions);
void print_energy(interactions_struct &interactions);

// Get the contact energy between 2 neighbouring sites with indices site1 and
// site2.
double get_contact_energy(state_struct& state,
                          int site1,
                          int site2,
                          interactions_struct& interactions,
                          geometry_space::Geometry& geometry);

// Get total energy of a given site, which is the sum of contact energies with
// its neighbours
double get_site_energy(state_struct& state,
                       interactions_struct& interactions,
                       geometry_space::Geometry& geometry,
                       int site_index);

// Get total energy of the system
double get_energy(state_struct& state,
                  interactions_struct& interactions,
                  geometry_space::Geometry& geometry);

std::ostream &operator<<(std::ostream &out, interactions_struct &interactions);
} // namespace lattice_particles_space

#endif
