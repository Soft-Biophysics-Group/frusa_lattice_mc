#include "utils.h"

namespace model_space{
  model_params::model_params(){
    /*
     * Populate the struct using the input JSON file
     */
  
    std::ifstream model_f;
    model_f.open("./model_params.json");
    if(!model_f){
      std::cerr << "Could not open JSON model parameters file" << std::endl;
      exit(1);
    }

    json json_model_params = json::parse(model_f);
  
    N          = json_model_params["N"].template get<int>();
    Np         = json_model_params["Np"].template get<int>();
    parameters = json_model_params["parameters"].template get<vec1d>();
  }

}
