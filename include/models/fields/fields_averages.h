#ifndef FIELDS_AVERAGES_HEADER_H
#define FIELDS_AVERAGES_HEADER_H

#include "fields_parameters.h"
#include "fields_state.h"
#include "fields_interactions.h"

#include <iostream>
#include <fstream>

#include "io_utils.h"
#include "vector_utils.h"

namespace fields_space{
  /*
   * Definitions required for the public routines of the model class
   */
  
  // Structure containing the containers for the average quantities:
  // e_av     - average MF energy 
  // e2_av    - average squared MF energy
  // s_av     - average MF entropy 
  // s2_av    - average squared MF entropy
  // f_av     - average MF free energy 
  // f2_av    - average squared MF free energy
  struct averages_struct{
    double e_av;
    double e2_av;
    double s_av;
    double s2_av;
    double f_av;
    double f2_av;
  };

  void initialize_averages(averages_struct &averages,
                           model_parameters_struct &parameters);
  
  void update_averages(averages_struct &averages, 
                       state_struct &state,
                       interactions_struct &interactions,
                       model_parameters_struct &parameters, 
                       double T);
  
  void save_averages(averages_struct &averages, 
                     state_struct &state,
                     model_parameters_struct &parameters, 
                     double T,int mcs_av);
  
  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */



}

#endif

