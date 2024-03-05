#include "model.h"

namespace model_space{
  model::model(){
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    initialize_state(state,parameters);
    
    //initialize_coupling_matrix(coupling_matrix);
    
    //get_energy(energy);

  }

  void model::print_model_state(){
    print_state(state, parameters);
  }

  void model::save_model_state(){
    save_state(state, parameters);
  }
}
