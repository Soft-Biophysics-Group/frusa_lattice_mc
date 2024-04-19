#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "vector_utils.h"
#include <array>
template <int M, int N>
using arr2i = std::array<std::array<int, N>,M>;

namespace geometry_space {
// Vectors of vectors giving the permutation indices for all the possible bond
// orientations.
struct bond_permutation_struct {
  static constexpr arr2i<6, 6> hexagonal{{{1, 2, 3, 4, 5, 6},
                                          {2, 3, 4, 5, 6, 1},
                                          {3, 4, 5, 6, 1, 2},
                                          {4, 5, 6, 1, 2, 3},
                                          {5, 6, 1, 2, 3, 4},
                                          {6, 1, 2, 3, 4, 5}}};
};

enum class lattices {
  one_dimension,
  square,
  triangular,
  cubic,
  bcc,
  fcc,
  n_lattices
};

class Geometry {
public:
  Geometry() = default;
  Geometry(lattices lattice, int lx, int ly, int lz);
  int get_neighbour(const int site_ind, const int bond_ind) const;
  int get_bond(const int site_1_ind, const int site_2_ind) const;
  int get_interaction_indices(const int site_1_ind, const int site_2_ind) const;
  double get_interaction(const int site_1_ind, const int site_2_ind,
                         const vec1d couplings) const;
private:
  lattices lattice_m {lattices::one_dimension};
  int lx_m {1};
  int ly_m {1};
  int lz_m {1};
  int n_neighbours {2};
};

} // namespace geometry_space
#endif
