#include "model.h"
#include "mc_routines.h"
#include <iostream>

int main(){

  std::cout << "Program runs!\n";
  
  model_space::model new_model;

  new_model.print_model_state();
  new_model.print_model_interactions();

  simulation_space::mc annealing;
  annealing.print_mc_parameters();
 
  annealing.t_scan(new_model);

  return 0;
}
