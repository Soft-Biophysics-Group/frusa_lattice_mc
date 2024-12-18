// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef DEFAULT_INTERACTIONS_HEADER_H
#define DEFAULT_INTERACTIONS_HEADER_H

#include "default_parameters.h"
#include "default_state.h"

#include <iostream>
#include <cmath>

#include "vector_utils.h"

namespace default_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing the characteristics of the model interactions:
  // delta  - energy gap
  // energy - current energy of the system
  struct interactions_struct{
    double delta;
    double energy;
  };

  // Calculate interaction characteristics of the current state of the system
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

  // Functions to calculate the energy of the system
  double get_energy(state_struct &state, double delta);
}
#endif
