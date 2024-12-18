// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef FIELDS_1D_HEADER_H
#define FIELDS_1D_HEADER_H

#include "utils.h"

namespace model_space{
  
  /*
   * Structure containing the fractional and total site concentrations:
   * psi_bar       - total particle density
   * concentration - fractional local concentrations of different orientations
   * local_density - local particle density
   */
  struct state_struct{
    double psi_bar;
    vec2d concentration;
    vec1d local_density;
  }

  /*
   * Structure containing information about the energetics of the system:
   * energy          - total energy of the system
   * coupling_matrix - two-site coupling map
   * T_model         - strength of the entropic mixing
   */
  struct interactions_struct{
    vec3d coupling_matrix;
    double T_model;
    double energy;
  }


  /*Definition of model class*/
  class model {
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

      /*Model parameter vector*/
      vec1d parameters;
      
      /*Interaction parameters*/
      double k11, k12, k21;

      /*Mean-field temperature (not annealing temperature!)*/
      double T_model;

      /*Vector containing the relative concentrations*/
      vec2d concentrations;

      /*Coupling matrix*/
      vec3d coupling_matrix;

      /*Energy of the system*/
      double energy;

      /*Threshold for smallest possible local density*/
      double eps=1e-8;
  
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

      /*Function to shift density from lattice site to another*/
      void shift_local_density(int,double);

      /*Calculate bounds for local density transfer*/
      void get_donor_bound(int, vec1d &, double &, int &);
      void get_acceptor_bound(int, vec1d &, double &, int &);

      /*Function to update the fractional concentrations on a given site*/
      void convert_concentrations(int,double);

      /*Calculate the change in the entropic contribution*/
      double get_entropy_shift(vec1d, int, double);
      double get_entropy_convert(vec1d, int, double);

      /*Function to calculate the density vector*/
      void update_psi(vec2d&);
    
    public:

      /*Class constructor*/
      model(const struct model_params &);
 
      /*Function to calculate the total energy of the system*/
      double get_energy();

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
