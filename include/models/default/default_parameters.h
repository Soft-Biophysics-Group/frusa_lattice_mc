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
  struct model_parameters_struct{
    model_parameters_struct();
    int N;
    double delta;
    EngineType rng;
    std::string initialize_option;
    std::string state_input;
  };
}  
#endif
