#include "geometry.h"
#include "triangular.h"
#include <stdexcept>

namespace geometry_space {
Geometry::Geometry(lattice_options lattice, int lx, int ly, int lz)
    : lattice_m{lattice}, lx_m{lx}, ly_m{ly}, lz_m{lz} {
  switch (lattice_m) {
  case lattice_options::triangular:
    n_neighbours_m = 6;
    n_orientations_m = 6;
    break;
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}

int Geometry::get_neighbour(const int site_ind, const int bond_ind) const {
  switch (lattice_m) {
  case lattice_options::triangular:
    return triangular_space::get_neighbour_triangular(site_ind, bond_ind, lx_m,
                                                      ly_m);
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}

int Geometry::get_bond(const int site_1_ind, const int site_2_ind) const {
  switch (lattice_m) {
  case lattice_options::triangular:
    return triangular_space::get_bond_triangular(site_1_ind, site_2_ind, lx_m,
                                                 ly_m);
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}

// TODO Uncomment if necessary, delete otherwise
//template <int N>
//arr1i<N>& Geometry::get_bond_permutation(const int site_1_ind, const int site_2_ind) const {
  //switch (lattice_m) {
  //case lattice_options::triangular:
    //return triangular_space::get_bond_permutation(site_1_ind, site_2_ind, lx_m,
                                                  //ly_m);
  //default:
    //throw(std::runtime_error("Invalid lattice option"));
  //}
//}

int Geometry::get_interaction_index(const int site_1_orientation,
                                    const int site_1_ind,
                                    const int site_2_orientation,
                                    const int site_2_ind) const {
  switch (lattice_m) {
  case lattice_options::triangular: {
    int bond{triangular_space::get_bond_triangular(site_1_ind, site_2_ind, lx_m,
                                                   ly_m)};
    return triangular_space::get_interaction_index(
        site_1_orientation, site_2_orientation, bond, n_orientations_m);
  }
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}

double Geometry::get_interaction(const int site_1_orientation,
                                 const int site_1_ind,
                                 const int site_2_orientation,
                                 const int site_2_ind,
                                 const vec1d flat_interaction_matrix) const {
  std::size_t interaction_index{
    static_cast<std::size_t>(
      get_interaction_index(
        site_1_orientation, site_1_ind, site_2_orientation, site_2_ind)
      )};
  return flat_interaction_matrix[interaction_index];
  }

int get_interaction_index(const int face_1, const int face_2,
                          const int n_orientations) {
  int interaction_index{};
  // Interaction atrix is encoded into 1D vector as upper-triangular
  // coefficients
  // Hijacking the hashing function Andrey already wrote
  if (face_1 <= face_2) {
    array_space::ij_to_r(interaction_index, face_1, face_2, n_orientations, 1);
  } else {
    array_space::ij_to_r(interaction_index, face_2, face_1, n_orientations, 1);
  }

  return interaction_index;
}
} // namespace geometry_space
