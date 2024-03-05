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


namespace model_space{
  

  /*
   * Generic routines that can be used by different models
   */

  /*Update average moments of energy*/
  void update_energy_moments(double,double,double);

  /*Update correlation function*/
  void update_correlation();

} 

#endif
