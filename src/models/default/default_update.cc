#include "default_update.h"

namespace default_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void update_system(state_struct &state, 
                     interactions_struct &interactions,
                     model_parameters_struct &parameters,
                     double T){
    
    int_dist particle_dist(0,state.N-1);
    real_dist uniform_dist(0,1);

    for(int r=0;r<state.N;r++){
      int index = particle_dist(parameters.rng);

      if(state.occupation[index]==1){
        state.occupation[index] = 0;
        state.average_occupation-= 1.0/state.N;
        interactions.energy-= interactions.delta;
      }
      else if(exp(-interactions.delta/T)>uniform_dist(parameters.rng)){
        state.occupation[index] = 1;
        state.average_occupation+= 1.0/state.N;
        interactions.energy+= interactions.delta;
      }
    } 
  }
  /* 
   * End of the required definitions for the model class
   */
}
