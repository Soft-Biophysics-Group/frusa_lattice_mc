// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#include "default_interactions.h"

namespace default_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_interactions(state_struct &state, 
                               interactions_struct &interactions,
                               model_parameters_struct &parameters){

    // Initialize model parameters from user input file
    interactions.delta = parameters.delta;

    // Calculate the energy of the initial state
    interactions.energy = get_energy(state, interactions.delta);
  }

  void print_interactions(state_struct &state,
                          interactions_struct &interactions){
    std::cout << "\n--------------------------------------------\n";
    std::cout << "Current interaction properties of the system\n";
    std::cout << "--------------------------------------------\n\n";
    
    std::cout << "The energy gap is\n\n";
    std::cout << "delta = " << interactions.delta << "\n\n";

    std::cout << "Current value of energy:\n\n";
    std::cout << "energy = " << interactions.energy/state.N << "\n\n";
  }

  void print_energy(state_struct &state, interactions_struct &interactions){
    std::cout << "\nE = " << interactions.energy/state.N << "\n\n";
  }
  /* 
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
  double get_energy(state_struct &state, double delta){

    double energy = 0;

    double e_levels[2] = {-delta/2,\
                           delta/2};

    for(int r=0;r<state.N;r++){
      energy+= e_levels[state.occupation[r]];
    }

    return energy;
  }
}
