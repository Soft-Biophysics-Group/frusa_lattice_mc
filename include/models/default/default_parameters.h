// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef DEFAULT_PARAMETERS_HEADER_H
#define DEFAULT_PARAMETERS_HEADER_H

#include <iostream>
#include <fstream>
#include <random>

#include <json.hpp>

#include "vector_utils.h"

typedef std::mt19937 EngineType;

using json = nlohmann::json;

namespace default_space{

  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing parameters used for model definition:
  // N                 - number of particles
  // delta             - energy gap between two states
  // rng               - random number generator
  // initialize_option - option string for choosing initialization function
  //                     current options are "from_file", "random", "uniform"
  // state_input       - if initialize_option is set to "from_file", this
  //                     string contains the location of the input structure
  // state_av_option   - boolean to determine if the simulation will store
  //                     MC averages of the occupation numbers
  // e_av_option       - boolean to determine if the simulation will store
  //                     MC averages of the first and second energy moments
  // state_av_output   - if state_av_option is true, gives the location for the
  //                     output file for MC averages of the occupation numbers
  // e_av_output       - if e_av_option is true, gives the location for the
  //                     output file for MC averages of the first and second
  //                     energy moments
  struct model_parameters_struct{
    model_parameters_struct();
    int N;
    double delta;
    EngineType rng;
    std::string initialize_option;
    std::string state_input;
    bool state_av_option;
    bool e_av_option;
    std::string state_av_output;
    std::string e_av_output;
  };
}
#endif
