#include "fields_averages.h"

namespace fields_space{
  
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
                       state_struct &state,
                       interactions_struct &interactions,
                       model_parameters_struct &parameters, 
                       double T){
  
    if(parameters.e_av_option==true){
      averages.e_av+= interactions.energy;
      averages.e2_av+= interactions.energy*interactions.energy;
      averages.s_av+= interactions.entropy;
      averages.s2_av+= interactions.entropy*interactions.entropy;
      averages.f_av+= interactions.free_energy;
      averages.f2_av+= interactions.free_energy*interactions.free_energy;
    }
  }

  void save_averages(averages_struct &averages, 
                     state_struct &state,
                     model_parameters_struct &parameters, 
                     double T, int mcs_av){
  
    if(parameters.e_av_option==true){
      averages.e_av/= mcs_av;
      averages.e2_av/= mcs_av;
      averages.s_av/= mcs_av;
      averages.s2_av/= mcs_av;
      averages.f_av/= mcs_av;
      averages.f2_av/= mcs_av;
      vec1d output_vec = {T,averages.e_av,averages.e2_av,
                            averages.s_av,averages.s2_av,
                            averages.f_av,averages.f2_av};
      std::string output_file = parameters.e_av_output\
                                +"esf_av_T_"+std::to_string(T)\
                                +".dat";
      io_space::save_vector(output_vec,7,output_file);
    }
  }

  /* 
   * End of the required definitions for the model class
   */
}

