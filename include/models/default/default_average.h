#ifndef DEFAULT_AVERAGE_HEADER_H
#define DEFAULT_AVERAGE_HEADER_H

#include "default_parameters.h"

#include <iostream>
#include <fstream>

#include "vector_utils.h"

namespace default_space{
  /*
   * Definitions required for the public routines of the model class
   */
  
  // Structure containing the containers for the average quantities:
  // average_state - average occupation number collected at fixed 
  //                 temperature
  struct averages_struct{
    double average_state;
  };

  void initialize_averages(averages_struct &averages, state_struct &state,
                           model_parameters_struct &parameters);
  
  void update_averages(averages_struct &averages, state_struct &state,
                       model_parameters_struct &parameters, double T);
  
  void save_averages(averages_struct &averages, state_struct &state,
                     model_parameters_struct &parameters, double T);
  
  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */



}

#endif
