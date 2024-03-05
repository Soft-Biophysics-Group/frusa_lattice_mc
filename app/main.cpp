#include "model.h"
//#include "mc_routines.h"
#include <iostream>

int main(){


  std::cout << "Program runs!\n";
  
  model_space::model_params parameters;

  std::cout << parameters.Lx << "\n";
  std::cout << parameters.Ly << "\n";
  std::cout << parameters.Lz << "\n";
  std::cout << parameters.Np << "\n";
  
  //simulation_space::mc_params parameters;
  //new_model.print_state();
  //new_model.print_energy();

  //simulation_space::mc<model_space::model> annealing(new_model, mc_params_1d);

  //annealing.t_scan();

  return 0;
}
