#include "model.h"

/*At compile time, select the model to be used in the simulation*/
//#if defined PARTICLES_1D_OPTION
//  #include "1d_particles.h"
//#if defined FIELDS_1D_OPTION
//  #include "1d_fields.h"
//#endif

namespace model_space{
  model::model(const model_params &parameters) :
    Lx(parameters.Lx),
    Ly(parameters.Ly),
    Lz(parameters.Lz),
    Np(parameters.Np),
    couplings(parameters.couplings),
    rng(parameters.rng)
  {
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    /*Calculate the total number of lattice sites*/
    N = Lx*Ly*Lz;

    //initialize_system(state);
    
    //initialize_coupling_matrix(coupling_matrix);
    
    //get_energy(energy);

  }
  
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
  
    Lx         = json_model_params["Lx"].template get<int>();
    Ly         = json_model_params["Ly"].template get<int>();
    Lz         = json_model_params["Lz"].template get<int>();
    Np         = json_model_params["Np"].template get<int>();
    couplings  = json_model_params["couplings"].template get<vec1d>();
  }
}
