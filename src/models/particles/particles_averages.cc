#include "particles_averages.h"
#include "io_utils.h"

namespace particles_space {

  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_averages(averages_struct &averages,
                           model_parameters_struct &parameters){

    if(parameters.e_av_option==true){
      averages.e_av = 0.0;
      averages.e2_av = 0.0;
    }
  }

  void update_averages(averages_struct &averages,
                       [[maybe_unused]] state_struct &state,
                       interactions_struct &interactions,
                       model_parameters_struct &parameters,
                       [[maybe_unused]] double T) {

    if(parameters.e_av_option==true){
      averages.e_av+= interactions.energy;
      averages.e2_av+= interactions.energy*interactions.energy;
    }
  }

  void save_averages(averages_struct &averages,
                     [[maybe_unused]] state_struct &state,
                     model_parameters_struct &parameters,
                     double T, int mcs_av){

    if(parameters.e_av_option==true){
      averages.e_av/= mcs_av;
      averages.e2_av/= mcs_av;
      vec1d output_vec = {T, averages.e_av, averages.e2_av};
      std::string output_file =
          std::to_string(static_cast<int>(parameters.e_av_output)) +
          "esf_av_T_" + std::to_string(T) + ".dat";
      io_space::save_vector(output_vec,7,output_file);
    }
  }

  /*
   * End of the required definitions for the model class
   */
} // lattice_particles_space
