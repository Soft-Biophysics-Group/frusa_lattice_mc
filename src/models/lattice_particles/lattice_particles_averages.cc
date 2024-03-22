#include "lattice_particles_averages.h"

namespace lattice_particles_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_averages(averages_struct &averages,
                           model_parameters_struct &parameters){
  
    if(parameters.state_av_option==true){
      averages.state_av = 0;
    }

    if(parameters.e_av_option==true){
      averages.e_av = 0;
      averages.e2_av = 0;
    }
  }

  void update_averages(averages_struct &averages, 
                       state_struct &state,
                       interactions_struct &interactions,
                       model_parameters_struct &parameters, 
                       double T){
  
    if(parameters.state_av_option==true){
      averages.state_av+= state.average_occupation;
    }

    if(parameters.e_av_option==true){
      averages.e_av+= interactions.energy;
      averages.e2_av+= interactions.energy*interactions.energy;
    }

  }

  void save_averages(averages_struct &averages, 
                     state_struct &state,
                     model_parameters_struct &parameters, 
                     double T, int mcs_av){
  
    if(parameters.state_av_option==true){
      averages.state_av/= mcs_av;
      vec1d output_vec = {T, averages.state_av};
      std::string output_file = parameters.state_av_output\
                                +"state_av_T_"+std::to_string(T)\
                                +".dat";
      io_space::save_vector(output_vec,2,output_file);
    }

    if(parameters.e_av_option==true){
      averages.e_av/= mcs_av;
      averages.e2_av/= mcs_av;
      vec1d output_vec = {T, averages.e_av,averages.e2_av};
      std::string output_file = parameters.state_av_output\
                                +"e_av_T_"+std::to_string(T)\
                                +".dat";
      io_space::save_vector(output_vec,3,output_file);
    }
  }

  /* 
   * End of the required definitions for the model class
   */
}
