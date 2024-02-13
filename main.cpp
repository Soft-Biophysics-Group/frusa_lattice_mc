#include "models.h"

int main(){

  simulation::model_data model_data_1d;

  model_data_1d.N   = 20;
  model_data_1d.Np  = 15;

  model_data_1d.k11 = -1;
  model_data_1d.k12 = 0;
  model_data_1d.k21 = 0;

  simulation::mc_data mc_data_1d;

  //mc_data_1d.T   = 0.01;
  
  static std::random_device dev;
  model_data_1d.rng = EngineType(dev());
  
  simulation::particles test(model_data_1d,mc_data_1d);

  test.print_state();
  //test.save_state("structure.dat","../Results/");
  test.print_energy();

  test.update_state(0.01);
  test.print_state();
  test.print_energy();

}
