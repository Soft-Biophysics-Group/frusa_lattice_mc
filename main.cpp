#include "models.h"

int main(){

  simulation::model_data model_data_1d;

  model_data_1d.N   = 20;
  model_data_1d.Np  = 5;

  model_data_1d.k11 = -1;
  model_data_1d.k12 = 0;
  model_data_1d.k21 = 0;

  simulation::particles test(model_data_1d);

  test.print_state();
  test.print_energy();

}
