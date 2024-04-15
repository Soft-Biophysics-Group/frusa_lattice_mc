// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#include "vector_utils.h"

namespace array_space {

int hash_3integers(int x1, int x2, int y, int x_range, int y_range) {
  return (x1 - 1) * x_range * y_range + (x2 - 1) * y_range + y;
}

void ij_to_r(int &r, int i, int j, int Lx) { r = i + Lx * j; }

void r_to_ij(int r, int &i, int &j, int Lx) {
  j = r / Lx;
  i = r - Lx * j;
}

void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly) {
  r = i + Lx * j + Lx * Ly * k;
}

void r_to_ijk(int r, int &i, int &j, int &k, int Lx, int Ly) {
  k = r / (Lx * Ly);
  j = (r - Lx * Ly * k) / Lx;
  i = r - Lx * j - Lx * Ly * k;
}

int mod(int a, int b) { return (a % b + b) % b; }
} // namespace array_space
