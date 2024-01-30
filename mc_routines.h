#ifndef MC_HEADER_H
#define MC_HEADER_H

#include "utils.h"

namespace simulation{
  class mc : public lattice_system::particles {
    private:
      /*Define random number distributions*/
      std::random_device dev;
      EngineType engine;
      int_dist particle_random;
      real_dist acceptance;

      /*
       * Simulation parameters
       */

      /*Number of MC steps for equilibration and averaging*/
      int mcs_eq, mcs_av;

      /*Real and log temperature*/
      double T, pT;

      /*Initial/final temperatures, step size*/
      double Ti, Tf, dT;

      /*Number of temperature steps*/
      int Nt;

      /*Average energy/energy^2, heat capacity*/
      double eav, e2av, cv;

      /*Average particle density*/
      vec2d psi_av;


      /*Options for data collections*/

      /*Collect average energy and heat capacity?*/
      bool ecv_option;

      /*Collect average particle density?*/
      bool psi_option;

      /*Save state configurations at intermediate temperatures?*/
      bool checkpoint;

    public:

      /*Class constructor*/
      mc();
  };

}

#endif
