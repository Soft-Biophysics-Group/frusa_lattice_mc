// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef IO_UTILS_HEADER_H
#define IO_UTILS_HEADER_H

#include "vector_utils.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace io_space{
  /*
   * Routines to assist with reading/writing data to/from files
   */

  // Read data from file at some location and write it to an array.
  // N, N1, N2 are dimensions of the target array
  // Currently using overloading, may be changed in the future to avoid
  // repeating definitions
  void read_vector(vec1i &array, int N, std::string& location);
  void read_vector(vec1d &array, int N, std::string& location);
  void read_vector(vec2i &array, int N1, int N2, std::string& location);
  void read_vector(vec2d &array, int N1, int N2, std::string& location);

  // Save data to a file at specified location
  void save_vector(vec1i &array, int N, std::string& location);
  void save_vector(vec1d &array, int N, std::string& location);
  void save_vector(vec2i &array, int N1, int N2, std::string& location);
  void save_vector(vec2d &array, int N1, int N2, std::string& location);
}

#endif
