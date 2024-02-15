#include "fields.h"
#include "mc_routines.h"

int main(){

  model_space::model_data model_data_1d;

  model_data_1d.N   = 20;
  model_data_1d.Np  = 5;

  model_data_1d.k11 = -1;
  model_data_1d.k12 = 0;
  model_data_1d.k21 = 0;

  simulation_space::mc_data mc_data_1d;

  mc_data_1d.mcs_eq = 1000;
  mc_data_1d.mcs_av = 1;
  mc_data_1d.Ti = 1;
  mc_data_1d.Tf = 0.01;
  mc_data_1d.Nt = 10;
  mc_data_1d.cooling_schedule = simulation_space::linear;
  mc_data_1d.checkpoint = true;
  mc_data_1d.checkpoint_address = "../Results/test/checkpoints/";

  static std::random_device dev;
  model_data_1d.rng = EngineType(dev());
  
  model_space::fields test(model_data_1d);

  test.print_state();
  test.print_energy();

  //test.update_state(0.01);
  
  simulation_space::mc<model_space::fields> annealing(test, mc_data_1d);

  annealing.t_scan();

  return 0;
}
