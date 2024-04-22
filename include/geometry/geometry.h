#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "vector_utils.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <array>
#include <string>

namespace geometry_space {
// Vectors of vectors giving the permutation indices for all the possible bond
// orientations.
enum lattice_options {
  one_dimension,
  square,
  triangular,
  cubic,
  bcc,
  fcc,
  n_lattices
};

static const inline std::array<std::string, lattice_options::n_lattices>
    lattice_str_arr{"chain", "square", "triangular", "cubic", "bcc", "fcc"};

lattice_options get_lattice_from_str(std::string& lattice_str);

class Geometry {
public:
  Geometry() = default;
  Geometry(lattice_options lattice, int lx, int ly, int lz=1);
  Geometry(std::string& geometry_input);
  int get_neighbour(const int site_ind, const int bond_ind) const;
  int get_bond(const int site_1_ind, const int site_2_ind) const;
  template <int N>
    arr1i<N>& get_bond_permutation(const int site_1_ind, const int site_2_ind) const;
  int get_interaction_index(const int site_1_orientation, const int site_1_ind,
                            const int site_2_orientation,
                            const int site_2_ind) const;
  double get_interaction(const int site_1_orientation, const int site_1_ind,
                         const int site_2_orientation, const int site_2_ind,
                         const vec1d flat_interaction_matrix) const;

private:
  lattice_options lattice_m {lattice_options::one_dimension};
  int lx_m {1};
  int ly_m {1};
  int lz_m {1};
  int n_neighbours_m {2};
  int n_orientations_m {2};
  void set_lattice_properties();
};

int get_interaction_index(const int face_1, const int face_2,
                          const int n_orientations);

} // namespace geometry_space
#endif
