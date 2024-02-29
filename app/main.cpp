#include "model.h"
#include "mc_routines.h"

int main(){

  model_space::model new_model = model_space::initialize_model();
  
  simulation_space::mc_params mc_params_1d;
  
  new_model.print_state();
  new_model.print_energy();

  simulation_space::mc<model_space::model> annealing(new_model, mc_params_1d);

  annealing.t_scan();

  return 0;
}
