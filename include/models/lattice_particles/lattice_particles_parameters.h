#ifndef LATTICEPARTICLES_PARAMETERS_HEADER_H
#define LATTICEPARTICLES_PARAMETERS_HEADER_H

#include <fstream>
#include <iostream>
#include <random>

#include <json.hpp>

#include "vector_utils.h"

typedef std::mt19937 EngineType;

using json = nlohmann::json;

namespace lattice_particles_space {
using ContactMap = vec1d;
enum mc_moves {
  swap_empty_full_enum,
  swap_full_full_enum,
  rotate_enum,
  mutate_enum,
  n_enum_moves
};
using move_probas_arr =
    std::array<double, static_cast<int>(mc_moves::n_enum_moves)>;


/*
 * Definitions required for the public routines of the model class
 */

// TODO Add try/catches to this to make sure that we define the right
// quantitites

// Structure containing parameters used for model definition:
// n_types           - number of different particle types
// n_orientations    - number of orientations a particle can take
// lx, ly, lz        - dimensions of the lattice
// n_particles       - number of particles
// couplings         - Flattened array of contact energies between each
//                     possible pair of faces
// rng               - random number generator
// initialize_option - option string for choosing initialization function
//                     current options are "from_file", "random"
// state_input       - if initialize_option is set to "from_file", this
//                     string contains the location of the input structure
struct model_parameters_struct {
  model_parameters_struct(const std::string &input_file);
  model_parameters_struct()
      : model_parameters_struct("./input/model_params.json"){};
  int n_types{};
  int n_orientations{};
  int lx{};
  int ly{};
  int lz{};
  vec1i n_particles{};
  ContactMap couplings{};
  EngineType rng{};
  std::string initialize_option{};
  std::string state_input{};
  move_probas_arr move_probas{};
  // TODO Add code to get option from json
  bool e_av_option {true};
  bool e_av_output {true};
};

std::ostream &operator<<(std::ostream &out, model_parameters_struct &params);

// Definition and choice of Monte Carlo move probabilities

// Returns a vector assigning the probability of each of the moves in
// mc_moves being picked at every MC step.
// mc_json_file is the same file used for in the engine header.
// The file should include a vector of n_enum_moves doubles summing to one,
// each of which corresponds to the probability of the associated move (in
// enum order) being picked.
move_probas_arr get_move_probas(const std::string &mc_json_file);

} // namespace lattice_particles_space
#endif
