#ifndef GEOMETRY_TRIANGULAR_H
#define GEOMETRY_TRIANGULAR_H

#include "vector_utils.h"
#include <unordered_map>
namespace geometry_space {
namespace triangular_space {
struct bond_permutation_struct {
  static constexpr arr2i<6, 6> bond_permutations{{{0, 1, 2, 3, 4, 5},
                                                  {5, 0, 1, 2, 3, 4},
                                                  {4, 5, 0, 1, 2, 3},
                                                  {3, 4, 5, 0, 1, 2},
                                                  {2, 3, 4, 5, 0, 1},
                                                  {1, 2, 3, 4, 5, 0}}};
  static constexpr arr2i<6,6> orientation_permutations = bond_permutations;
  static inline const std::unordered_map<vec1i, int> bond_indices{
      {{1, 0}, 0},  {{1, 1}, 1},   {{-1, 1}, 2},
      {{-1, 0}, 3}, {{-1, -1}, 4}, {{1, -1}, 5}};
};

int get_neighbour_triangular(const int site_ind, const int bond_ind,
                             const int lx, const int ly);

// Returns 6 if the sites are not neighbours
int get_bond_triangular(const int site_1_ind, const int site_2_ind,
                        const int lx, const int ly);

int get_face(const int orientation, const int bond);

int get_interaction_index(const int site_1_orientation,
                          const int site_2_orientation, const int bond_1,
                          const int n_orientations);

int get_opposite_face(const int face);
} // namespace triangular_space

} // namespace geometry_space

#endif
