// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef IO_VECTOR_HEADER_H
#define IO_VECTOR_HEADER_H

#include <vector>

typedef std::vector<int> vec1i;
typedef std::vector<std::size_t> vec1s;
typedef std::vector<double> vec1d;
using vec1b = std::vector<bool>;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

namespace array_space{
  /*
   * Routines to assist with manipulation of flattened arrays.
   * Flattened arrays are preferred over nested vector because during the
   * MC updates we only have to generate a single random index, instead of 3.
   */
  // Next few functions deal with going from vectors of parameters to hashed
  // integers for efficient contact energy retrieval

// Get a hashed contact integer from pair characteristics
// Hash 3 integers, the 2 first of which can take the same values, into a
// single number
int hash_3integers(int x1, int x2, int y, int x_range, int y_range);

// Next functions are for going back and forth between lattice coordinates
// and a single hashed site index

// Calculates the flat array index (r) given the 2D coordinates (i,j) and
// the dimensions of the system (Lx,Ly)
void ij_to_r(int &r, int i, int j, int Lx);

// Calculates the 2D coordinates (i,j) given the flat array inde (r) and
// the dimensions of the system
void r_to_ij(int r, int &i, int &j, int Lx);

// Calculates the flat array index (r) given the 3D coordinates (i,j,k) and
// the dimensions of the system (Lx,Ly,Lz)
void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly, int Lz);

// Calculates the 3D coordinates (i,j,k) given the flat array inde (r) and
// the dimensions of the system
void r_to_ijk(int r, int &i, int &j, int &k, int Lx, int Ly, int Lz);

// Calculates modulo of an integer a (positive or negative) with respect to
// another positive integer b
int mod(int a, int b);

template <typename T>
void print_vector(std::ostream &out, const std::vector<T> &vec) {
  for (T elem : vec)
    out << elem << ',';
  }

template <typename T, int N>
void print_array(std::ostream &out, const std::array<T, N> &vec) {
  for (T elem : vec)
    out << elem << ',';
}

} // namespace array_space


#endif
