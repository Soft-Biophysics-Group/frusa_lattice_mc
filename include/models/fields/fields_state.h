#ifndef FIELDS_STATE_HEADER_H
#define FIELDS_STATE_HEADER_H

#include "utils.h"

namespace fields_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing the fractional and total site concentrations:
  // ns            - number of orientation species
  // Lx, Ly, Lz    - dimensions of the lattice 
  // N             - number of lattice sites
  // rho_bar       - total particle density
  // concentration - fractional local concentrations of different orientations
  // local_density - local particle density
  struct state_struct{
    int ns;
    int Lx;
    int Ly;
    int Lz;
    int N;
    double rho_bar;
    vec2d concentration;
    vec1d local_density;
  }

  // Initialize the structural properties of the system, depending on the type 
  // of parameters.initialize_option
  void initialize_state(state_struct &state, 
                        model_parameters_struct parameters);

 
  // Print the current values of the structural properties of the system
  void print_state(state_struct state, model_parameters_struct parameters);

  // Save the fractional concentrations to a file 
  void save_state(state_struct state, model_parameters_struct parameters);
  
  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */

  // Various methods for state initialization
  void initialize_state_from_file(state_struct &state, 
                                  model_parameters_struct parameters);

  void initialize_state_random(state_struct &state, 
                               model_parameters_struct parameters);

  void initialize_state_uniform(state_struct &state);

  // Returns a vector with the positions of the neighbours of a specified site
  vec1i get_neighbours(int r, state_struct state);
}
