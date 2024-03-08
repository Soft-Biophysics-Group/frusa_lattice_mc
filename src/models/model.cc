#include "model.h"

namespace model_space{
  model::model(){
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    initialize_state(state,parameters);
    
    initialize_interactions(state,interactions,parameters);
    
  }

  void model::print_model_state(){
    print_state(state);
  }

  void model::save_model_state(){
    save_state(state, parameters);
  }
  
  void model::print_model_interactions(){
    print_interactions(state,interactions);
  }

}
