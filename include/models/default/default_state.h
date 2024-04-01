// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef DEFAULT_STATE_HEADER_H
#define DEFAULT_STATE_HEADER_H

#include "default_parameters.h"

#include <iostream>
#include <random>

#include "io_utils.h"
#include "vector_utils.h"

typedef std::mt19937 EngineType;
typedef std::uniform_int_distribution<int> int_dist;
typedef std::uniform_real_distribution<double> real_dist;

namespace default_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing the characteristics of the state of the system:
  // N                  - number of particles
  // average_occupation - average occupation number of the current state 
  // occupation         - array of state occupations
  struct state_struct{
    int N;
    double average_occupation;
    vec1i occupation;
  };

  // Initialize the structural properties of the system, depending on the type 
  // of parameters.initialize_option
  void initialize_state(state_struct &state, 
                        model_parameters_struct &parameters);

  // Print the state properties of the system
  void print_state(state_struct &state);

  // Save the occupation numbers to a file "state_output"
  void save_state(state_struct &state, std::string state_output);
  
  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */

  // Various methods for state initialization
  void initialize_state_from_file(state_struct &state, 
                                  model_parameters_struct &parameters);

  void initialize_state_random(state_struct &state, 
                               model_parameters_struct &parameters);

  void initialize_state_uniform(state_struct &state);
}
#endif
