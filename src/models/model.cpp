#include "model.h"

/*At compile time, select the model to be used in the simulation*/
#if defined PARTICLES_1D_OPTION
  #include "1d_particles.h"
#endif
//#ifdef FIELDS_1D_OPTION
//  #include "1d_fields.h"
//#endif


namespace model_space{
  model initialize_model(){
    /*
     * Function used to initialize the selected model
     */

    /*Initialize and fill model parameters structure using the input file*/
    model_params parameters;

    static std::random_device dev;
    parameters.rng = EngineType(dev());

    /*Initialize simulation model*/
    model m(parameters);

    return m;
  }
}
