#ifndef MODEL_HEADER_H
#define MODEL_HEADER_H

#include <iostream>
#include <fstream>

#include <random>
typedef std::mt19937 EngineType;
typedef std::uniform_int_distribution<int> int_dist;
typedef std::uniform_real_distribution<double> real_dist;

#include <vector>
typedef std::vector<int> vec1i;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

#include <json.hpp>
using json = nlohmann::json;


namespace model_space{
  /*Definition of model class*/
  class model {
    private:
      /*
       * Private variables
       */

      /*Dimensions of the lattice*/
      int Lx, Ly, Lz;

      /*Total number of lattice sites*/
      int N;

      /*Number of particles*/
      int Np;

      /*Vector describing the state of the system*/
      //state_vector state;

      /*Model parameter vector*/
      vec1d couplings;
      
      /*Coupling matrix*/
      vec3d coupling_matrix;

      /*Energy of the system*/
      double energy;

      /*Random number engine to be used by the class*/  
      EngineType rng;

      /*
       * Private routines of the class
       */

      /*Function to initialize the system of anisotropic particles*/
      void initialize();
    
    public:

      /*Class constructor*/
      model(const struct model_params &);
 
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
 
  struct model_params{
  /*Structure containing model parameters*/
    model_params();
    int N;
    int Np;
    vec1d couplings;
    EngineType rng;
  };

}
#endif
