#include "fields_parameters.h"

namespace fields_space{

  model_parameters_struct::model_parameters_struct(){
  
    std::ifstream model_input_f;
    model_input_f.open("./input/model_params.json");
    if(!model_input_f){
      std::cerr << "Could not open JSON model parameters file" << std::endl;
      exit(1);
    }

    json json_model_params = json::parse(model_input_f);
  
    ns         = json_model_params["ns"].template get<int>();
    Lx         = json_model_params["Lx"].template get<int>();
    Ly         = json_model_params["Ly"].template get<int>();
    Lz         = json_model_params["Lz"].template get<int>();
    Np         = json_model_params["Np"].template get<int>();
    couplings         = json_model_params["couplings"].template get<vec1d>();
    initialize_option = \
      json_model_params["initialize_option"].template get<std::string>();
    
    if(initialize_option=="from_file"){
      state_input = \
      json_model_params["state_input"].template get<std::string>();
    }
    
    std::random_device dev;
    EngineType engine(dev());

    rng = engine;
  }
}
