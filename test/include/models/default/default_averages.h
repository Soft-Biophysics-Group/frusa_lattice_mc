// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef DEFAULT_AVERAGES_HEADER_H
#define DEFAULT_AVERAGES_HEADER_H

#include "default_parameters.h"
#include "default_state.h"
#include "default_interactions.h"

#include <iostream>
#include <fstream>

#include "io_utils.h"
#include "vector_utils.h"

namespace default_space{
  /*
   * Definitions required for the public routines of the model class
   */
  
  // Structure containing the containers for the average quantities:
  // state_av - average occupation number collected at fixed 
  //                 temperature
  // e_av     - average energy 
  // e2_av    - average squared energy (used to calculate heat capacity)
  struct averages_struct{
    double state_av;
    double e_av;
    double e2_av;
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
