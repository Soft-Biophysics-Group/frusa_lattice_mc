#ifndef LATTICEPARTICLES_UPDATE_HEADER_H
#define LATTICEPARTICLES_UPDATE_HEADER_H

#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"
#include "lattice_particles_geometry.h"
#include "lattice_particles_interactions.h"
#include "json.hpp"

namespace lattice_particles_space{
  using json = nlohmann::json;
  /*
   * Definitions required for the public routines of the model class
   */

  // Function used to perform state.N updates of the system (considered a
  // single MC step)
  // state        - configuration of the system before the update
  // interactions - energetics of the system before the update
  // parameters   - used to access random number generator (parameters.rng)
  // T            - annealing temperature (not the same as T_model!)
  template <int N>
  void update_system(state_struct &state,
                     interactions_struct<N> &interactions,
                     model_parameters_struct &parameters,
                     double T);

  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */

  enum class mc_moves{
    swap_empty_full,
    swap_full_full,
    rotate,
    mutate,
    n_enum_moves
  };
  using move_probas_vec = std::array<double, static_cast<int>(mc_moves::n_enum_moves)>;

  // Returns a vector assigning the probability of each of the moves in
  // mc_moves being picked at every MC step.
  // mc_json_file is the same file used for in the engine header.
  // The file should include a vector of n_enum_moves doubles summing to one,
  // each of which corresponds to the probability of the associated move (in
  // enum order) being picked.
  move_probas_vec get_move_probas(std::string& mc_json_file);

  // Returns the index of a random lattice site containing a particle
  int select_random_full(state_struct& state, model_parameters_struct &parameters);

  // Same for a random lattice site containing no particle
  int select_random_empty(state_struct& state, model_parameters_struct &parameters);

  /*
   * IMPORTANT NOTE: The following functions only perform the candidate MC move
   * but do not calculate the resulting energy difference.
   */
  void swap_empty_full(state_struct &state, model_parameters_struct &parameters);
  void swap_full_full(state_struct &state, model_parameters_struct &parameters);
  void rotate(state_struct &state, model_parameters_struct &parameters);
  void mutate(state_struct &state, model_parameters_struct &parameters);
}

#endif
