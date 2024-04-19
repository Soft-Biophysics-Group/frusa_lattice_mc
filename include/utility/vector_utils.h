// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef IO_VECTOR_HEADER_H
#define IO_VECTOR_HEADER_H

#include <vector>
#include <array>

typedef std::vector<int> vec1i;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

template <int M, int N>
using arr2i = std::array<std::array<int, N>,M>;

namespace array_space{
  /*
   * Routines to assist with manipulation of flattened arrays.
   * Flattened arrays are preferred over nested vector because during the 
   * MC updates we only have to generate a single random index, instead of 3.
   */

  // Calculates the flat array index (r) given the 2D coordinates (i,j) and
  // the dimensions of the system (Lx,Ly)
  void ij_to_r(int &r, int i, int j, int Lx, int Ly);
  
  // Calculates the 2D coordinates (i,j) given the flat array inde (r) and
  // the dimensions of the system
  void r_to_ij(int r, int &i, int &j, int Lx, int Ly);

  // Calculates the flat array index (r) given the 3D coordinates (i,j,k) and
  // the dimensions of the system (Lx,Ly,Lz)
  void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly, int Lz);
  
  // Calculates the 3D coordinates (i,j,k) given the flat array inde (r) and
  // the dimensions of the system
  void r_to_ijk(int r, int &i, int &j, int &k, int Lx, int Ly, int Lz);

  // Calculates modulo of an integer a (positive or negative) with respect to
  // another positive integer b
  int mod(int a, int b);
}

#endif
