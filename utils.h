#ifndef MC_UTILS_HEADER_H
#define MC_UTILS_HEADER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>

typedef std::mt19937 EngineType;
typedef std::vector<int> vec1i;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

namespace lattice_system{
  class particles {
    private:
      
      /*Help randomize the seed*/
      std::random_device dev;

      /*
       * Private routines of the class
       */

      /*Function to initialize the system of anisotropic particles*/
      int initialize(int,int);
    
    public:
      
      /*State vector containing the positions and 
       * orientations of the particles */
      vec2i state;

      /*Coupling matrix*/
      vec3d coupling_matrix;

      /*Energy of the system*/
      double E;

      /*
       * Public routines of the class
       */

      /*Class constructor*/
      particles(int, int, double, double, double);

      /*Function to extract the states of neighbours of a given particle*/
      vec1i get_neighbours(int, int, int);

      /*Function to calculate the total energy of the system*/
      double get_energy(int, int);

      /*Function to calculate the Fourier Transform of the state vector*/
      int get_state_q();

  };
} 

#endif
