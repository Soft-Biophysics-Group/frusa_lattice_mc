// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "model.h"
#include "mc_routines.h"
#include <iostream>


// Temporary includes for testing
#include "vector_utils.h"

int main(){
  std::cout << "Program runs!\n";

  //particles_space::move_probas_arr test{
      //particles_space::get_move_probas("input/model_params.json")};
  //array_space::print_array(std::cout, test);

  model_space::model new_model;

  new_model.print_model_state();
  new_model.print_model_interactions();

  simulation_space::mc annealing;
  annealing.print_mc_parameters();

  annealing.t_scan(new_model);

  return 0;
}
