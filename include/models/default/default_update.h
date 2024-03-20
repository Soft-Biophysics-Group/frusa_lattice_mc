#ifndef DEFAULT_UPDATE_HEADER_H
#define DEFAULT_UPDATE_HEADER_H

#include "default_parameters.h"
#include "default_state.h"
#include "default_interactions.h"

namespace default_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Function used to perform state.N updates of the system (considered a
  // single MC step)
  // state        - configuration of the system before the update
  // interactions - energetics of the system before the update
  // parameters   - used to access random number generator (parameters.rng)
  // T            - annealing temperature (not the same as T_model!)
  void update_system(state_struct &state, 
                     interactions_struct &interactions,
                     model_parameters_struct &parameters,
                     double T);

  /*
   * End of the required definitions for the model class
   */
}

#endif
