// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef FCC_GEOMETRY_HEADER_H
#define FCC_GEOMETRY_HEADER_H

#include "vector_utils.h"
#include <iostream>
#include <map>

using BondIndexMap = std::map<std::array<int,3> , unsigned int>;
namespace geometry_space {
  namespace fcc_space {

    struct bond_structure {
      static inline const vec2i bond_permutation = {
        // We group bonds/faces into triangles, which do not share vertices,
        // and define permutations using the 8 C_3 [111] rotations and
        // 3 C_2 [100] rotations.
        {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11 }, // E
        {  2,  0,  1,  8,  6,  7, 11,  9, 10,  5,  3,  4 }, // C_3(1)
        { 11,  9, 10,  5,  3,  4,  2,  0,  1,  8,  6,  7 }, // C_3(2)
        {  5,  3,  4, 11,  9, 10,  8,  6,  7,  2,  0,  1 }, // C_3(3)
        {  8,  6,  7,  2,  0,  1,  5,  3,  4, 11,  9, 10 }, // C_3(4)
        {  1,  2,  0, 10, 11,  9,  4,  5,  3,  7,  8,  6 }, // C_3^2(1)
        {  7,  8,  6,  4,  5,  3, 10, 11,  9,  1,  2,  0 }, // C_3^2(2)
        { 10, 11,  9,  1,  2,  0,  7,  8,  6,  4,  5,  3 }, // C_3^2(3)
        {  4,  5,  3,  7,  8,  6,  1,  2,  0, 10, 11,  9 }, // C_3^2(4)
        {  3,  4,  5,  0,  1,  2,  9, 10, 11,  6,  7,  8 }, // C_2(1)
        {  9, 10, 11,  6,  7,  8,  3,  4,  5,  0,  1,  2 }, // C_2(2)
        {  6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5 }  // C_2(3)
      }; // bond_permutation

      static inline const vec2i bond_array {
        {  1,  0,  0 },
        {  0,  1,  0 },
        {  0,  0,  1 },
        { -1,  0,  0 },
        { -1,  0,  1 },
        { -1,  1,  0 },
        {  0,  1, -1 },
        {  1,  0, -1 },
        {  0,  0, -1 },
        {  0, -1,  1 },
        {  0, -1,  0 },
        {  1, -1,  0 },
      }; // bond_array

      static inline const BondIndexMap bond_index {
        { {  1,  0,  0 }, 0 },
        { {  0,  1,  0 }, 1 },
        { {  0,  0,  1 }, 2 },
        { { -1,  0,  0 }, 3 },
        { { -1,  0,  1 }, 4 },
        { { -1,  1,  0 }, 5 },
        { {  0,  1, -1 }, 6 },
        { {  1,  0, -1 }, 7 },
        { {  0,  0, -1 }, 8 },
        { {  0, -1,  1 }, 9 },
        { {  0, -1,  0 }, 10 },
        { {  1, -1,  0 }, 11 }
      }; // bond_index
    }; // bond_structure

    void get_bond();
    void get_neighbour();
    void get_interaction();
    void get_interaction_indices();

  } // fcc_space
} // geometry_space

#endif // FCC_GEOMETRY_HEADER_H
