#ifndef MODELS_HEADER_H
#define MODELS_HEADER_H

#include "utils.h"

namespace simulation{
  class particles {
    private:
      
      /*Help randomize the seed*/
      std::random_device dev;
      
      /*
       * Private variables
       */

      /*Lattice size*/
      int N;

      /*Number of particles*/
      int Np;

      /*Interaction parameters*/
      double k11, k12, k21;

      /*State vector containing the positions and orientations of the 
        particles*/
      vec2i state;

      /*Coupling matrix*/
      vec3d coupling_matrix;

      /*Energy of the system*/
      double energy;

      /*
       * Private routines of the class
       */

      /*Function to initialize the system of anisotropic particles*/
      void initialize();
      
      /*Function to extract the states of neighbours of a given particle*/
      vec1i get_neighbours(int);

      /*Function to calculate the total energy of the system*/
      double get_energy();

      /*Function to calculate the density vector*/
      void update_psi(vec2d&);
    
    public:

      /*Class constructor*/
      particles(const struct model_data &);
 
      /*
       * Required public routines of the class
       */

      /*Print the current state of the system*/
      void print_state();

      /*Print the current energy of the system*/
      void print_energy();
       
  };
} 

#endif
