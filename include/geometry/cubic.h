#ifndef GEOMETRY_CUBIC_H
#define GEOMETRY_CUBIC_H

#include "vector_utils.h"
#include <map>
using BondIndexMap = std::map<std::array<int,3> , int>;

namespace geometry_space {
namespace cubic_space {
struct bond_struct {
  // See notes from 25/04, to be put cleaner
  static inline const vec2i bond_permutation{{
    {0, 1, 2, 3, 4, 5}, //E
    {1, 6, 2, 3, 0, 4}, // C4(z)-
    {2, 1, 5, 0, 4, 3}, // C4(y)+
    {3, 1, 0, 5, 4, 3}, // C4(y)-
    {4, 0, 2, 3, 5, 1}, // C4(x)+
    {5, 1, 3, 2, 4, 0}, // C4(y)+^2
  }}; //bond_permutation

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
};
static constexpr int n_neighbours{6};
static constexpr int n_orientations{6};
} // namespace cubic_space
} // namespace geometry_space

#endif
