#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "vector_utils.h"
#include <array>

namespace geometry_space {
// Vectors of vectors giving the permutation indices for all the possible bond
// orientations.
enum class lattice_options {
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
  Geometry(lattice_options lattice, int lx, int ly, int lz=1);
  int get_neighbour(const int site_ind, const int bond_ind) const;
  int get_bond(const int site_1_ind, const int site_2_ind) const;
  int get_interaction_indices(const int site_1_ind, const int site_2_ind) const;
  double get_interaction(const int site_1_ind, const int site_2_ind,
                         const vec1d couplings) const;
private:
  lattice_options lattice_m {lattice_options::one_dimension};
  int lx_m {1};
  int ly_m {1};
  int lz_m {1};
  int n_neighbours {2};
};

} // namespace geometry_space
#endif
