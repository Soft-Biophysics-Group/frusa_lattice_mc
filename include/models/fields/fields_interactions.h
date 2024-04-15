// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef FIELDS_INTERACTIONS_HEADER_H
#define FIELDS_INTERACTIONS_HEADER_H

#include "fields_parameters.h"
#include "fields_state.h"
#include "fields_geometry.h"

#include <iostream>
#include <cmath>

#include "vector_utils.h"

namespace fields_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing the characteristics of the model interactions:
  // coupling_matrix - interaction map between different fields
  // T_model         - strength of the mean-field entropic interactions
  // energy          - current energy of the system
  // entropy         - current mean-field entropy of the system
  // free_energy     - energy + entropy
  struct interactions_struct{
    vec3d coupling_matrix;
    double T_model;
    double energy;
    double entropy;
    double free_energy;
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
  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */

  // Functions to calculate mean-field energy and entropy of the system
  // eps=1e-8 is used to determine threshold for calculating logs
  double get_energy(state_struct &state, vec3d coupling_matrix);
  double get_entropy(state_struct &state, double T_model, double eps = 1e-8);

  // Function to update energy after concentration transfer from donor to the
  // acceptor
  // field.
  // r_d             - lattice position of the donor
  // r_a             - lattice position of the acceptor
  // index_d         - donor field component
  // index_a         - acceptor field component
  // dc              - amount of concentration transferred
  double get_energy_change(int r_d, int r_a, int index_d, int index_a,
                           double dc, state_struct &state,
                           interactions_struct &interactions);

  // Function to update entropy after concentration is increased or decreased
  // on a single site (local density not conserved)
  // r             - lattice position of the updated concentration field
  // index         - updated field component
  // dc            - amound of concentration gained or lost
  // eps           - threshold for log calculation
  double get_entropy_change_shift(int r, int index, double dc,
                                  state_struct &state,
                                  interactions_struct &interactions,
                                  double eps=1e-8);

  // Function to update entropy after concentration is transferred between
  // components of the field on a single site (local density conserved)
  // See arguments for get_entropy_change_shift and
  // index_d - component of the donor field
  // index_a - component of the acceptor field
  double get_entropy_change_convert(int r, int index_d, int index_a, double dc,
                                    state_struct &state,
                                    interactions_struct &interactions,
                                    double eps=1e-8);

  // Function to update mean-field interaction properties
  // dE - change in energy
  // dS - change in entropy
  // dF - change in free energy
  void update_interactions(double dE, double dS, double dF,
                           interactions_struct &interactions);
}
#endif
