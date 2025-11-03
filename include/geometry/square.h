// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef GEOMETRY_SQUARE_H
#define GEOMETRY_SQUARE_H

#include "vector_utils.h"
#include <map>

using BondIndexMap = std::map<std::array<int,3> , int>;

namespace geometry_space {
namespace square_space {
struct bond_struct {
  static inline const vec2i bond_permutation{{
    { 0, 1, 2, 3 },
    { 1, 2, 3, 0 },
    { 2, 3, 0, 1 },
    { 3, 0, 1, 2 }
  }}; // bond_permutation

  static inline const vec2i bond_array{
    {  1,  0, 0 },
    {  0,  1, 0 },
    { -1,  0, 0 },
    {  0, -1, 0 }
  }; // bond_array

  static inline const BondIndexMap bond_index{
    { { 1,  0, 0}, 0 },
    { { 0,  1, 0}, 1 },
    { {-1,  0, 0}, 2 },
    { { 0,  -1, 0}, 3 },
  }; //bond_index

  static inline const vec1i opposite_bonds {2, 3, 0, 1};
}; // bond_struct
static constexpr int n_neighbours {4};
static constexpr int n_orientations {4};
  } // square_space
} // geometry_space
 
#endif //SQUARE_GEOMETRY_HEADER_H
