#include "model.h"
//#include "mc_routines.h"
#include <iostream>

int main(){

  std::cout << "Program runs!\n";
  
  model_space::model new_model;

  new_model.print_model_state();
  new_model.print_model_interactions();

  //simulation_space::mc<model_space::model> annealing(new_model, mc_params_1d);

  //annealing.t_scan();

  return 0;
}
