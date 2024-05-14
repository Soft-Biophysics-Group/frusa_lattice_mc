// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef CHAIN_GEOMETRY_HEADER_H
#define CHAIN_GEOMETRY_HEADER_H

#include "vector_utils.h"
#include <map>
using BondIndexMap = std::map<std::array<int,3> , int>;

namespace geometry_space {
  namespace chain_space {

    struct bond_struct{
      static inline const vec2i bond_permutation {{
          { 0, 1 }, { 1, 0 }
        }};
        //constexpr unsigned int bond_permutation[2][2] = {
            //{0, 1}, {1, 0}}; // bond_permutation

        static inline const vec2i bond_array = {{1, 0, 0},
                                                {-1, 0, 0}}; // bond_array
        // constexpr int bond_array[2][3] = {{1, 0, 0}, {-1, 0, 0}}; //
        // bond_array

        static inline const BondIndexMap bond_index{
            {{ 1, 0, 0}, 0},
            {{-1, 0, 0}, 1},
        };
        //std::map<int[3], unsigned int> bond_index{
            //{{1, 0, 0}, 0},
            //{{-1, 0, 0}, 1},
        //}; // bond_index
      }; // bond_structure

    //void get_bond();
    //void get_neighbour();
    //void get_interaction();
    //void get_interaction_indices();
    static constexpr int n_neighbours{2};
    static constexpr int n_orientations{2};

  } // chain_space
} // geometry_space

#endif //CHAIN_GEOMETRY_HEADER_H
