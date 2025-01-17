// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include <iostream>
#include <CLI11.hpp>

#include "mc_routines.h"
#include "model.h"

int main(int argc, char** argv)
{
  std::cout << "Program runs!\n";

  // Define the command-line arguments we can parse
  CLI::App app{"Performs a Monte-Carlo simulated annealing on a lattice."};
  argv = app.ensure_utf8(argv);

  std::string model_params_file = "input/model_params.json";
  app.add_option("-m,--model-params",
                 model_params_file,
                 "JSON input file for model parameters");

  std::string mc_params_file = "input/mc_params.json";
  app.add_option("-M,--mc-params",
                 mc_params_file,
                 "JSON input file for Monte-Carlo annealing parameters");

  CLI11_PARSE(app, argc, argv);

  // Create model
  model_space::model new_model {model_params_file};

  new_model.print_model_state();
  new_model.print_model_interactions();

  // Create MC engine and perform annealing
  simulation_space::mc annealing;
  annealing.print_mc_parameters();

  annealing.t_scan(new_model);

  return 0;
}
