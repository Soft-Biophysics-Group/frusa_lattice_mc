#ifndef LATTICEPARTICLES_AVERAGES_H
#define LATTICEPARTICLES_AVERAGES_H

#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"
#include "lattice_particles_interactions.h"

#include <iostream>
#include <fstream>

namespace lattice_particles_space {
  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing the containers for the average quantities:
  // e_av     - average energy
  // e2_av    - average squared energy
  struct averages_struct{
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

}

#endif
