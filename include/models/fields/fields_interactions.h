#ifndef FIELDS_INTERACTIONS_HEADER_H
#define FIELDS_INTERACTIONS_HEADER_H

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
  // coupling_matrix - interaction map between different particles
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
  
  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
  
  // Calculate mean-field energy and entropy of the system
  // eps=1e-8 is used to determine threshold for calculating logs
  double get_energy(state_struct &state, vec3d coupling_matrix);
  double get_entropy(state_struct &state, double T_model, double eps=1e-8);
  
  // Calculate the difference in energy and entropy after local updates
  void get_energy_diff();
  void get_entropy_diff();


}
#endif
