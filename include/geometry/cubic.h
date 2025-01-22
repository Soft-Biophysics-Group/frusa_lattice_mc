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
  static inline const vec2i bond_permutation {{
    {
      0,  1,  2,  3,
      4,  5,  6,  7,
      8,  9, 10, 11,
      12, 13, 14, 15,
      16, 17, 18, 19,
      20, 21, 22, 23
    }, //E
    {
      4,  5,  6,  7,
      12, 15, 14, 13,
      9, 10, 11,  8,
      16, 17, 18, 19,
      0,  3,  2,  1,
      21, 22, 23, 20
    }, // C4Z-
    {
       8,  9, 10, 11,
       7,  4,  5,  6,
      14, 13, 12, 15,
      20, 21, 22, 23,
      19, 16, 17, 18,
       2,  1,  0,  3
    }, // C4Y+
    {
      12, 15, 14, 13,
      16, 19, 18, 17,
      10, 11,  8,  9,
       0,  3,  2,  1,
       4,  7,  6,  5,
      22, 23, 20, 21
    }, // C2Z
    {
      16, 19, 18, 17,
       0,  1,  2,  3,
      11, 8,  9,  10,
       4,  7,  6,  5,
      12, 13, 14, 15,
      23, 20, 21, 22
    }, // C4Z+
    {
      20, 23, 22, 21,
      19, 18, 17, 16,
      12, 15, 14, 13,
       8, 11, 10,  9,
       7,  6,  5,  4,
       0,  3,  2,  1
    }, // C2Z * C4Y+
  }};  // bond_permutation

  //static inline const std::unordered_map<vec1i, int> bond_indices{
  static inline const vec2i bond_array{
    { 1,  0,  0},
    { 0,  1,  0},
    { 0,  0,  1},
    {-1,  0,  0},
    { 0, -1,  0},
    { 0,  0, -1},
  }; //bond_directions

  static inline const BondIndexMap bond_index{
    { { 1,  0,  0}, 0 },
    { { 0,  1,  0}, 1 },
    { { 0,  0,  1}, 2 },
    { {-1,  0,  0}, 3 },
    { { 0, -1,  0}, 4 },
    { { 0,  0, -1}, 5 },
  }; //bond_index

  static inline const vec1i opposite_bonds {3, 4, 5, 0, 1, 2};
};  // bond_structure
static constexpr int n_neighbours {6};
static constexpr int n_orientations {24};
}  // namespace cubic_space
}  // namespace geometry_space

#endif
