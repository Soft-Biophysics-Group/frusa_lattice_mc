#ifndef FIELDS_UPDATE_HEADER_H
#define FIELDS_UPDATE_HEADER_H

#include "fields_state.h"
#include "fields_interactions.h"

namespace fields_space{
  /*
   * Definitions required for the public routines of the model class
   */
 void update_state(state_struct &state, 
                   interactions_struct &interactions,
                   model_parameters_struct parameters,
                   double T);

  /*
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
 
 // Function to shift density from one lattice site to another
 void shift_local_density(state_struct &state,
                          interactions_struct &interactions,
                          model_parameters_struct &parameters,
                          double T);

 /*Function to update the fractional concentrations on a given site*/
 void convert_concentrations(int,double);
  
}

#endif
