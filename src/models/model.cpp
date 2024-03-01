#include "model.h"

/*At compile time, select the model to be used in the simulation*/
#if defined PARTICLES_1D_OPTION
  #include "1d_particles.h"
#elif defined FIELDS_1D_OPTION
  #include "1d_fields.h"
#endif

namespace model_space{
  model::model(const model_params &parameters) :
    N(parameters.N),
    Np(parameters.Np),
    couplings(parameters.couplings),
    rng(parameters.rng)
  {
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    initialize_system(state);
    
    initialize_coupling_matrix(coupling_matrix);
    
    get_energy(energy);

  }
}
