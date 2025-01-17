// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include <iostream>

#include "mc_routines.h"
#include "model.h"

int main()
{
  std::cout << "Program runs!\n";

  // particles_space::move_probas_arr test{
  // particles_space::get_move_probas("input/model_params.json")};
  // array_space::print_array(std::cout, test);

  model_space::model new_model;
  geometry_space::Geometry geometry {
      geometry_space::lattice_options::triangular, 40, 40, 1};
  int test_neighbour {geometry.get_neighbour(0, 0)};
  std::cout << "Test neighbour is" << test_neighbour << "\n";

  new_model.print_model_state();
  new_model.print_model_interactions();

  simulation_space::mc annealing;
  annealing.print_mc_parameters();

  annealing.t_scan(new_model);

  return 0;
}
