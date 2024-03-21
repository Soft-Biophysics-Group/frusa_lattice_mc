#include "default_average.h"

namespace default_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_averages(averages_struct &averages, state_struct &state,
                           model_parameters_struct &parameters){
  
    if(parameters.average_state_option==true){
      averages.average_state = 0;
    }
  }

  void update_averages(averages_struct &averages, state_struct &state,
                           model_parameters_struct &parameters, double T){
  
    if(parameters.average_state_option==true){
      averages.average_state+= state.average_occupation;
    }
  }

  void save_averages(averages_struct &averages, state_struct &state,
                     model_parameters_struct &parameters, 
                     double T, int mcs_av){
  
    if(parameters.average_state_option==true){
      averages.average_state/= mcs_av;
     
    }
  }

  /* 
   * End of the required definitions for the model class
   */
}
