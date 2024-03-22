#include "default_parameters.h"

namespace default_space{

  model_parameters_struct::model_parameters_struct(){
  
    std::ifstream model_input_f;
    model_input_f.open("./input/model_params.json");
    if(!model_input_f){
      std::cerr << "Could not open JSON model parameters file" << std::endl;
      exit(1);
    }

    json json_model_params = json::parse(model_input_f);
  
    N     = json_model_params["N"].template get<int>();
    delta = json_model_params["delta"].template get<double>();
    initialize_option = \
      json_model_params["initialize_option"].template get<std::string>();
    
    if(initialize_option=="from_file"){
      state_input = \
      json_model_params["state_input"].template get<std::string>();
    }

    state_av_option = \
      json_model_params["state_av_option"].template get<bool>();
    e_av_option     = json_model_params["e_av_option"].template get<bool>();
     
    if(state_av_option==true){
      state_av_output = \
      json_model_params["state_av_output"].template get<std::string>();
    }
    if(e_av_option==true){
      e_av_output = \
      json_model_params["e_av_output"].template get<std::string>();
    }

    std::random_device dev;
    EngineType engine(dev());

    rng = engine;
  }
}
