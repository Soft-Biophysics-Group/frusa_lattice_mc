#ifndef LATTICEPARTICLES_INTERACTIONS_H
#define LATTICEPARTICLES_INTERACTIONS_H

#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"
#include "lattice_particles_geometry.h"

#include <iostream>
#include <cmath>

#include "vector_utils.h"

namespace lattice_particles_space {
  // Structure containing the characteristics of the model interactions:
  // energy          - current energy of the system
  // free_energy     - energy + entropy
  struct interactions_struct{
    ContactMap coupling_matrix;
    double energy;
  };

  // Calculate interactions characteristics of the current state of the system
  void initialize_interactions(state_struct &state,
                               interactions_struct &interactions,
                               model_parameters_struct &parameters);

  // Print the summary of the interactions characteristics
  void print_interactions(state_struct &state,
                          interactions_struct &interactions);

  void print_energy(state_struct &state,
                    interactions_struct &interactions);

  double get_site_energy(int site_index, state_struct &state,
                         ContactMap coupling_matrix);

  // Function to update energy after concentration transfer from donor to
  // the acceptor field. r_d             - lattice position of the donor r_a
  // - lattice position of the acceptor index_d         - donor field
  // component index_a         - acceptor field component dc              -
  // amount of concentration transferred
  double get_energy_change(int r_d, int r_a, int index_d, int index_a,
                           double dc, state_struct &state,
                           interactions_struct &interactions);

} // lattice_particles_space

#endif
