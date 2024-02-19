#ifndef UTILS_HEADER_H
#define UTILS_HEADER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <time.h>
#include <iomanip>

#include <json.hpp>

using json = nlohmann::json;

typedef std::mt19937 EngineType;
typedef std::uniform_int_distribution<int> int_dist;
typedef std::uniform_real_distribution<double> real_dist;

typedef std::vector<int> vec1i;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

namespace model_space{
  /*
   * Data structures used to store the relevant parameters
   */
  
  /*Model parameters*/
  struct model_params{
    model_params();
    int N;
    int Np;
    vec1d parameters;
    EngineType rng;
  };

  /*
   * Generic routines that can be used by different models
   */

  /*Update average moments of energy*/
  void update_energy_moments(double,double,double);

  /*Update correlation function*/
  void update_correlation();

} 

#endif
