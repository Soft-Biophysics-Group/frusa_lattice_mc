// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef FIELDS_UPDATE_HEADER_H
#define FIELDS_UPDATE_HEADER_H

#include "fields_parameters.h"
#include "fields_state.h"
#include "fields_geometry.h"
#include "fields_interactions.h"

namespace fields_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Function used to perform state.N updates of the system (considered a
  // single MC step)
  // state        - configuration of the system before the update
  // interactions - energetics of the system before the update
  // parameters   - used to access random number generator (parameters.rng)
  // T            - annealing temperature (not the same as T_model!)
  void update_system(state_struct &state, 
                     interactions_struct &interactions,
                     model_parameters_struct &parameters,
                     double T);

  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
 
  // Function to shift density from one lattice site to another
  // See arguments of update_state
  // eps - threshold for density transfer (transfer everything below eps)
  void shift_local_density(state_struct &state,
                           interactions_struct &interactions,
                           model_parameters_struct &parameters,
                           double T, double eps=1e-8);

  // Function to update the fractional concentrations on a given site
  // See arguments for update_state and shift_local_density
  void convert_concentrations(state_struct &state,
                              interactions_struct &interactions,
                              model_parameters_struct &parameters,
                              double T, double eps=1e-8);  

  // Function to randomly select a field component at a given lattice site
  // r          - lattice position of the field
  // bound      - field component must be !=bound in order to be 
  //              a selection candidate
  // state      - current state of the system
  // parameters - used to extract random number generator (parameters.rng) 
  int select_element(int r, double bound, 
                     state_struct &state, 
                     model_parameters_struct &parameters);
}

#endif
