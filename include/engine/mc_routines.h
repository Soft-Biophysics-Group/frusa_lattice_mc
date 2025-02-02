// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef MC_HEADER_H
#define MC_HEADER_H

#include <iostream>
#include <cmath>

#include <json.hpp>

#include "vector_utils.h"
#include "model.h"

using json = nlohmann::json;

namespace simulation_space{

  /*Monte Carlo parameters*/
  struct mc_parameters_struct{
    mc_parameters_struct(std::string& mc_input);
    int mcs_eq {};
    int mcs_av {};
    double Ti {};
    double Tf {};
    int Nt {};
    std::string cooling_schedule {};
    bool checkpoint_option {};
    std::string checkpoint_address {};
    std::string final_structure_address {};
  };

  class mc {
    private:
      /*
       * Private variables
       */

      // Simulation parameters
      mc_parameters_struct parameters;

      // Integer option for annealing schedule
      int cooling_option {};

      //Array with annealing temperatures
      vec1d T_array {};

    public:

      // Class constructor
      mc(std::string& mc_input);

      // Prints user-defined MC parameters
      void print_mc_parameters();

      // MC annealing
      void extracted();
      void t_scan(model_space::model &simulation_system);

      // MC simmulation at a fixed temperature T
      void mc_simulate(model_space::model &simulation_model, double T);
  };
}

#endif
