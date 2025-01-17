#ifndef GEOMETRY_TRIANGULAR_H
#define GEOMETRY_TRIANGULAR_H

#include "vector_utils.h"
#include <map>
using BondIndexMap = std::map<std::array<int,3> , int>;

namespace geometry_space {
namespace triangular_space {
struct bond_struct {
  static inline const vec2i bond_permutation{{
      {0, 1, 2, 3, 4, 5},
      {1, 2, 3, 4, 5, 0},
      {2, 3, 4, 5, 0, 1},
      {3, 4, 5, 0, 1, 2},
      {4, 5, 0, 1, 2, 3},
      {5, 0, 1, 2, 3, 4},
  }}; //bond_permutation

  //static inline const std::unordered_map<vec1i, int> bond_indices{
  static inline const vec2i bond_array{
    { 1,  0, 0},
    { 0,  1, 0},
    {-1,  1, 0},
    {-1,  0, 0},
    { 0, -1, 0},
    { 1, -1, 0},
  }; //bond_directions

  static inline const BondIndexMap bond_index{
    { { 1,  0, 0}, 0 },
    { { 0,  1, 0}, 1 },
    { {-1,  1, 0}, 2 },
    { {-1,  1, 0}, 3 },
    { { 0, -1, 0}, 4 },
    { { 1, -1, 0}, 5 },
  }; //bond_index

  static inline const vec1i opposite_bonds {3, 4, 5, 0, 1, 2};
};  // bond_structure
static constexpr int n_neighbours {6};
static constexpr int n_orientations {6};
}  // namespace triangular_space
}  // namespace geometry_space

#endif
