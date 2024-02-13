#ifndef MODELS_HEADER_H
#define MODELS_HEADER_H

#include "utils.h"

namespace simulation{
  class particles {
    private:
      /*
       * Private variables
       */

      /*Lattice size*/
      int N;

      /*Number of particles*/
      int Np;

      /*Interaction parameters*/
      double k11, k12, k21;

      /*Vectors containing the positions and orientations of the particles*/
      vec1i positions;
      vec1i orientations;

      /*Coupling matrix*/
      vec3d coupling_matrix;

      /*Energy of the system*/
      double energy;
  
      /*
       * Pseudorandom number generator definitions
       */

      /*Define random number engine to be used by the class*/  
      EngineType rng;

      /*Define random distributions used by the class*/
      real_dist uniform_dist;
      int_dist  particle_dist;
      int_dist  empty_dist;
      int_dist  binary_dist;

      /*
       * Private routines of the class
       */

      /*Function to initialize the system of anisotropic particles*/
      void initialize();
      
      /*Function to extract the states of neighbours of a given particle*/
      vec1i get_neighbours(int);

      /*Function to calculate the total energy of the system*/
      double get_energy();

      /*Function to update the position of the chosen particle*/
      void update_position(int,double);

      /*Function to update the orientation of the chosen particle*/
      void update_orientation(int,double);

      /*Function to calculate the density vector*/
      void update_psi(vec2d&);
    
    public:

      /*Class constructor*/
      particles(const struct model_data &, const struct mc_data &);
 
      /*
       * Required public routines of the class
       */

      /*Print the current state of the system*/
      void print_state();

      /*Save the current state of the system to a file*/
      void save_state(std::string, std::string);

      /*Print the current energy of the system*/
      void print_energy();

      /*Update the state of the system*/
      void update_state(double);

      /*Update the selected simulation averages*/
      void update_averages();

      /*Save the selected simulation averages to a file*/
      void save_averages();
       
  };
} 

#endif
