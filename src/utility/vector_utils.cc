// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "vector_utils.h"

namespace array_space {

void ij_to_r(int& r, int i, int j, int Lx, [[maybe_unused]] int Ly)
{
  r = i + Lx * j;
}

void r_to_ij(int r, int &i, int &j, int Lx, [[maybe_unused]] int Ly) {
  j = r / Lx;
  i = r - Lx * j;
}

void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly,
              [[maybe_unused]] int Lz) {
  r = i + Lx * j + Lx * Ly * k;
}

void r_to_ijk(const int r,
              int& i,
              int& j,
              int& k,
              const int Lx,
              const int Ly,
              [[maybe_unused]] const int Lz)
{
  k = r / (Lx * Ly);
  j = (r - Lx * Ly * k) / Lx;
  i = r - Lx * j - Lx * Ly * k;
}

int mod(int a, int b) { return (a % b + b) % b; }

//template <typename T>
//void print_vector(std::ostream &out, const std::vector<T> &vec) {
  //out << '{';
  //for (int coeff : vec) {
    //out << coeff << ", ";
  //}
  //out << "}";
//}

//template <typename T, std::size_t N>
//void print_array(std::ostream &out, const std::array<T, N> &vec) {
  //out << '{';
  //for (int coeff : vec) {
    //out << coeff << ", ";
  //}
  //out << "}";
//}

//template <typename T>
//void print_vector_2d(std::ostream &out, const std::vector<std::vector<T>> &vec_vec) {
  //out << "{\n";
  //for (vec1i& vec : vec_vec) {
    //print_vector<T>(out, vec);
    //out << ",\n";
  //}
  //out << "}";
//}
} // namespace array_space
