// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "vector_utils.h"

namespace array_space{

  void ij_to_r(int &r, int i, int j, int Lx, int Ly){
    r = i + Lx*j;
  }

  void r_to_ij(int r, int &i, int &j, int Lx, int Ly){
    j = r/Lx;
    i = r-Lx*j;
  }

  void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly, int Lz){
    r = i + Lx*j + Lx*Ly*k;
  }

  void r_to_ijk(int r, int &i, int &j, int &k, int Lx, int Ly, int Lz){
    k = r/(Lx*Ly);
    j = (r-Lx*Ly*k)/Lx;
    i = r-Lx*j-Lx*Ly*k;
  }

  int mod(int a, int b){
    return (a % b + b) % b;
  }
}
