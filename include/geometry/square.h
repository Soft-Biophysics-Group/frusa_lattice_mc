// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef SQUARE_GEOMETRY_HEADER_H
#define SQUARE_GEOMETRY_HEADER_H

#include <iostream>

namespace geometry_space {
  namespace square_space {
  
    struct bond_structure {
      bond_structure();
      
      constexpr unsigned int bond_permutation[4][4] = {
        { 0, 1, 2, 3 },
        { 1, 2, 3, 0 },
        { 2, 3, 0, 1 },
        { 3, 0, 1, 2 }
      }; // bond_permutation

      constexpr int bond_array[4][3] = {
        {  1,  0, 0 },
        {  0,  1, 0 },
        { -1,  0, 0 },
        {  0, -1, 0 }
      }; // bond_array

      std::map<int[3], unsigned int> bond_index {
        { 0, {  1,  0, 0 } },
        { 1, {  0,  1, 0 } },
        { 2, { -1,  0, 0 } },
        { 3, {  0, -1, 0 } }
      }; // bond_index
    }; // bond_structure

    void get_bond();
    void get_neighbour();
    void get_interaction();
    void get_interaction_indices();

  } // square_space
} // geometry_space
 
#endif //SQUARE_GEOMETRY_HEADER_H
