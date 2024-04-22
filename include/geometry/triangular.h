#ifndef GEOMETRY_TRIANGULAR_H
#define GEOMETRY_TRIANGULAR_H

#include "vector_utils.h"
#include <map>
namespace geometry_space {
namespace triangular_space {
struct bond_permutation_struct {
  static inline vec2i bond_permutations{{{0, 1, 2, 3, 4, 5},
                                         {1, 2, 3, 4, 5, 0},
                                         {2, 3, 4, 5, 0, 1},
                                         {3, 4, 5, 0, 1, 2},
                                         {4, 5, 0, 1, 2, 3},
                                         {5, 0, 1, 2, 3, 4}}};
  //static inline const std::unordered_map<vec1i, int> bond_indices{
  static inline vec2i bond_directions{{1, 0, 0},  {1, 1, 0},   {-1, 1, 0},
                                      {-1, 0, 0}, {-1, -1, 0}, {1, -1, 0}};
};
} // namespace triangular_space
} // namespace geometry_space

#endif
