// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef FCC_GEOMETRY_HEADER_H
#define FCC_GEOMETRY_HEADER_H

#include <iostream>

namespace geometry_space {
  namespace fcc_space {
  
    struct bond_structure {
      constexpr unsigned int bond_permutation[12][12] = \
      {{1,2,3,4,5,6,7,8,9,10,11,12},
       {},
       {},
       {},
       {},
       {},
       {},
       {},
       {},
       {},
       {},
       {}}; // bond_permutation

      constexpr std::unordered_map<unsigned int, int[3]> bond_index {
        {0,},
        {1,},
        {2,},
        {3,},
        {4,},
        {5,},
        {6,},
        {7,},
        {8,},
        {9,},
        {10,},
        {11,}
      }; // bond_index
    
    }; // bond_structure

    void get_bond();
    void get_neighbour();
    void get_interaction();
    void get_interaction_indices();

  } // fcc_space
} // geometry_space
 
#endif // FCC_GEOMETRY_HEADER_H
