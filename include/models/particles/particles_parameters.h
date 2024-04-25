#ifndef PARTICLES_PARAMETERS_HEADER_H
#define PARTICLES_PARAMETERS_HEADER_H

#include <fstream>
#include <iostream>
#include <random>

#include <json.hpp>

#include "vector_utils.h"

typedef std::mt19937 EngineType;

using json = nlohmann::json;

namespace particles_space {

// Specific type alias for the contact map/couplings, in case we want to change
// it down the line
// TODO Remove this type

// Enum for the different types of possible particle moves.
enum mc_moves {
  swap_empty_full,
  swap_full_full,
  rotate,
  mutate,
  rotate_and_move,
  n_enum_moves
};

static const inline std::array<std::string, n_enum_moves> mc_moves_str{
    "swap_empty_full",
    "swap_full_full",
    "rotate",
    "mutate",
    "rotate_and_move"
};

// User-supplied array of probabilities of selecting each type of move during
// lattice update
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
// n_particles       - number of particles of each type
// couplings         - Flattened array of contact energies between each
//                     possible pair of faces
// rng               - random number generator
// initialize_option - option string for choosing initialization function
//                     current options are "from_file", "random"
// state_input       - if initialize_option is set to "from_file", this
//                     string contains the location of the input structure
// move_probas       - user-supplied array of various moves' probabilities of
//                     being picked during update
struct model_parameters_struct {
  model_parameters_struct(const std::string &model_input_file);
  model_parameters_struct()
      : model_parameters_struct("./input/model_params.json"){};
  int n_types{};
  vec1i n_particles{};
  vec1d couplings{};
  EngineType rng{};
  std::string initialize_option{};
  std::string state_input{};
  move_probas_arr move_probas{};
  // TODO Add code to get option from json; understand what these do
  bool e_av_option {true};
  bool e_av_output {true};
};

std::ostream &operator<<(std::ostream &out, model_parameters_struct &params);

// Definition and choice of Monte Carlo move probabilities

// Returns a vector assigning the probability of each of the moves in
// mc_moves being picked at every MC step.
// The json file should include a vector of n_enum_moves doubles summing to one,
// each of which corresponds to the probability of the associated move (in
// enum order) being picked.
move_probas_arr get_move_probas(const std::string &model_input_file);

} // namespace lattice_particles_space
#endif
