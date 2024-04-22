// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef CHAIN_GEOMETRY_HEADER_H
#define CHAIN_GEOMETRY_HEADER_H

#include <iostream>

namespace geometry_space {
  namespace chain_space {
  
    struct bond_structure {
      bond_structure();
      
      constexpr unsigned int bond_permutation[2][2] = {
        { 0, 1 },
        { 1, 0 }
      }; // bond_permutation

      constexpr int bond_array[2][3] = {
        {  1, 0, 0 },
        { -1, 0, 0 }
      }; // bond_array

      std::map<int[3], unsigned int> bond_index {
        { {  1, 0, 0 }, 0 },
        { { -1, 0, 0 }, 1 },
      }; // bond_index
    }; // bond_structure

    void get_bond();
    void get_neighbour();
    void get_interaction();
    void get_interaction_indices();

  } // chain_space
} // geometry_space
 
#endif //CHAIN_GEOMETRY_HEADER_H
