#ifndef GEOMETRY_TRIANGULAR_H
#define GEOMETRY_TRIANGULAR_H

#include "vector_utils.h"
#include <unordered_map>
namespace geometry_space {
struct bond_permutation_struct {
  static constexpr arr2i<6, 6> bond_permutations{{{1, 2, 3, 4, 5, 6},
                                                  {2, 3, 4, 5, 6, 1},
                                                  {3, 4, 5, 6, 1, 2},
                                                  {4, 5, 6, 1, 2, 3},
                                                  {5, 6, 1, 2, 3, 4},
                                                  {6, 1, 2, 3, 4, 5}}};
  static inline const std::unordered_map<vec1i, int> bond_indices{
      {{1, 0}, 0},  {{1, 1}, 1},   {{-1, 1}, 2},
      {{-1, 0}, 3}, {{-1, -1}, 4}, {{1, -1}, 5}};
};

namespace triangular_space {
int get_neighbour_triangular(const int site_ind, const int bond_ind,
                             const int lx, const int ly);
} // namespace triangular_space

} // namespace geometry_space

#endif
