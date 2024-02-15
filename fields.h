#ifndef FIELDS_HEADER_H
#define FIELDS_HEADER_H

#include "utils.h"

namespace model_space{
  class fields {
    private:
      /*
       * Private variables
       */

      /*Lattice size*/
      int N;

      /*Number of particles*/
      int Np;

      /*Total density of particles*/
      double psi_bar;

      /*Interaction parameters*/
      double k11, k12, k21;

      /*Vector containing the relative concentrations*/
      vec2d concentrations;

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
      int_dist  binary_dist;

      /*
       * Private routines of the class
       */

      /*Function to initialize the system of anisotropic fields*/
      void initialize();
      
      /*Function to extract the states of neighbours of a given site*/
      vec2i get_neighbours(int);

      /*Function to calculate the total energy of the system*/
      double get_energy();

      /*Function to shift density from lattice site to another*/
      void shift_local_density(int,double);

      /*Calculate bounds for local density transfer*/
      void get_donor_bound(int, double &, int &);
      void get_acceptor_bound(int, double &, int &);

      /*Function to update the fractional concentrations on a given site*/
      void update_local_concentrations(int,double);

      /*Function to calculate the density vector*/
      void update_psi(vec2d&);
    
    public:

      /*Class constructor*/
      fields(const struct model_data &);
 
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

      /*Initialize the containers to store the selected averages*/
      void initialize_averages();

      /*Update the selected simulation averages*/
      void update_averages(double);

      /*Save the selected simulation averages to the corresponding files*/
      void save_averages();
       
  };
} 

#endif
