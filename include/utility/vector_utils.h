// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef IO_VECTOR_HEADER_H
#define IO_VECTOR_HEADER_H

#include <array>
#include <vector>
#include <iostream>

typedef std::vector<int> vec1i;
typedef std::vector<std::size_t> vec1s;
typedef std::vector<double> vec1d;
typedef std::vector<std::size_t> vec1s;
typedef std::vector<std::vector<int>> vec2i;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

template<std::size_t N>
using arr1i = std::array<int, N>;
template<std::size_t M, std::size_t N>
using arr2i = std::array<arr1i<N>, M>;

namespace array_space {
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
void r_to_ijk(
    const int r, int& i, int& j, int& k, const int Lx, const int Ly, int Lz);

// Calculates modulo of an integer a (positive or negative) with respect to
// another positive integer b
int mod(int a, int b);

template <typename T>
void print_vector(std::ostream &out, const std::vector<T> &vec) {
  out << '{';
  for (T coeff : vec) {
    out << coeff << ", ";
  }
  out << "}";
}


template <typename T>
void print_vector_2d(std::ostream &out, const std::vector<std::vector<T>> &vec_vec) {
  out << "{\n";
  for (const vec1i& vec : vec_vec) {
    print_vector<T>(out, vec);
    out << ",\n";
  }
  out << "}";
}

template <typename T, std::size_t N>
void print_array(std::ostream &out, const std::array<T, N> &vec) {
  out << '{';
  for (T coeff : vec) {
    out << coeff << ", ";
  }
  out << "}";
}


} // namespace array_space

#endif
