#ifndef GEOMETRY_CUBIC_H
#define GEOMETRY_CUBIC_H

#include "vector_utils.h"
#include <map>
using BondIndexMap = std::map<std::array<int,3> , int>;

/**
 * When considering three-dimensional particles, we have to keep in mind that
 * the interactions for two particles in contact does not only depend on which
 * faces touch each other.
 * Instead, the relative particle orientation has to also be taken into
 * account. This can easily be understood by playing with model cubes!
 * We then have to take into account the 24 possible orientations cubes can
 * have, and how these permute when looked at through different bond
 * directions.
 */

namespace geometry_space {
namespace cubic_space {
struct bond_struct {
  static inline const vec2i bond_permutation{{
      {
            0,  1,  2,  3,
            4,  5,  6,  7,
            8,  9,  10, 11,
            12, 13, 14, 15,
            16, 17, 18, 19,
            20, 21, 22, 23,
      }, // E
      {
            4,   5,  6,  7,
            23, 20, 21, 22,
            11, 10,  9,  8,
            15, 12, 13, 14,
            0,  1,  2,  3,
            18, 19, 16, 17,
      }, // C4(z)+
      {
             8,  9, 10, 11,
             7,  4,  5,  6,
            21, 22, 23, 20,
             0,  1,  2,  3,
            17, 18, 19, 16,
            12, 13, 14, 15,
      }, // C4(y)+
      {
            12, 13, 14, 15, 
             5,  6,  7,  4,
             0,  1,  2,  3,
            21, 22, 23, 20,
            19, 16, 17, 18,
             8,  9, 10, 11,
      }, // C4(y)-
      {
            16, 17, 18, 19,
             0,  1,  2,  3,
            11,  8,  9, 10,
            13, 14, 15, 12,
            23, 20, 21, 22,
             6,  7,  4,  5,
      }, // C4(z)+
      {
            21, 22, 23, 20,
             6,  7,  4,  5,
            12, 13, 14, 15,
            8,  9,  10, 11,
            18, 19, 16, 17,
            0,  1,  2,  3,
      }, // C4(y)+^2
  }}; // bond_permutation

  //static inline const std::unordered_map<vec1i, int> bond_indices{
  static inline const vec2i bond_array{
    { 1,  0,  0},
    { 0,  1,  0},
    { 0,  0,  1},
    { 0,  0, -1},
    { 0, -1,  0},
    {-1,  0,  0},
  }; //bond_directions

  static inline const BondIndexMap bond_index{
    { { 1,  0,  0}, 0 },
    { { 0,  1,  0}, 1 },
    { { 0,  0,  1}, 2 },
    { { 0,  0, -1}, 3 },
    { { 0, -1,  0}, 4 },
    { {-1,  0,  0}, 5 },
  }; //bond_index
  
  static inline const vec1i opposite_bonds {5, 4, 3, 2, 1, 0};
};  // bond_structure
static constexpr int n_neighbours {6};
static constexpr int n_orientations {24};
}  // namespace cubic_space
}  // namespace geometry_space

#endif
