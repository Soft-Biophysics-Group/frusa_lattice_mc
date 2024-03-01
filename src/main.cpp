#include "particles.h"
#include "mc_routines.h"

int main(){

  model_space::model_params model_params_1d;
  
  static std::random_device dev;
  model_params_1d.rng = EngineType(dev());
  
  model_space::particles test(model_params_1d);

  simulation_space::mc_params mc_params_1d;
  
  test.print_state();
  test.print_energy();

  simulation_space::mc<model_space::particles> annealing(test, mc_params_1d);

  annealing.t_scan();

  return 0;
}
