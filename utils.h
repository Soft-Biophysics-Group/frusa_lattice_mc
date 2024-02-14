#ifndef UTILS_HEADER_H
#define UTILS_HEADER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <time.h>

typedef std::mt19937 EngineType;
typedef std::uniform_int_distribution<int> int_dist;
typedef std::uniform_real_distribution<double> real_dist;

typedef std::vector<int> vec1i;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

namespace simulation{
  /*
   * Data structures used to store the relevant simulation parameters
   */
  
  /*Model parameters*/
  struct model_data{
    int N;
    int Np;
    double k11;
    double k12;
    double k21;
    EngineType rng;
  };

  /*Options for the cooling schedule*/
  enum cooling_option {exponential, linear};

  /*Monte Carlo parameters*/
  struct mc_data{
    int mcs_eq;
    int mcs_av;
    double Ti;
    double Tf;
    int Nt;
    cooling_option cooling_schedule;
    bool checkpoint;
    std::string checkpoint_address;
  };
     
} 

#endif
