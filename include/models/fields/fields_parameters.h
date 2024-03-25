#ifndef FIELDS_PARAMETERS_HEADER_H
#define FIELDS_PARAMETERS_HEADER_H

#include <iostream>
#include <fstream>
#include <random>

#include <json.hpp>

#include "vector_utils.h"

typedef std::mt19937 EngineType;

using json = nlohmann::json;

namespace fields_space{

  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing parameters used for model definition:
  // ns                - number of orientation species
  // Lx, Ly, Lz        - dimensions of the lattice 
  // Np                - number of particles
  // couplings         - array containing the coupling matrix on +x bond
  // rng               - random number generator
  // initialize_option - option string for choosing initialization function
  //                     current options are "from_file", "random", "uniform"
  // state_input       - if initialize_option is set to "from_file", this 
  //                     string contains the location of the input structure
  // e_av_option       - boolean to determine if the simulation will store 
  //                     MC averages of the first and second MF E/S/F moments
  // e_av_output       - if e_av_option is true, gives the location for the  
  //                     output file for MC averages of the first and second 
  //                     MF E/S/F moments moments
  struct model_parameters_struct{
    model_parameters_struct();
    int ns;
    int Lx;
    int Ly;
    int Lz;
    int Np;
    vec1d couplings;
    EngineType rng;
    std::string initialize_option;
    std::string state_input;
    bool e_av_option;
    std::string e_av_output;
  };
}  
#endif
