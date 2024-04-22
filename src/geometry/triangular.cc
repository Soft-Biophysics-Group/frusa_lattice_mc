#include "triangular.h"
#include "vector_utils.h"


namespace geometry_space {

// Forward declaration to avoid circular dependency
int get_interaction_index(const int face_1, const int face_2,
                          const int n_orientations);

namespace triangular_space {
int get_neighbour_triangular(const int site_ind, const int bond_ind,
                             const int lx, const int ly) {

  int i{};
  int j{};
  int neigh_ind{0};

  array_space::r_to_ij(site_ind, i, j, lx, ly);

  int im = array_space::mod(i - 1, lx);
  int ip = array_space::mod(i + 1, lx);
  int jm = array_space::mod(j - 1, ly);
  int jp = array_space::mod(j + 1, ly);

  switch (bond_ind) {
  case 0:
    array_space::ij_to_r(neigh_ind, ip, j, lx, ly);
    break;
  case 1:
    array_space::ij_to_r(neigh_ind, ip, jp, lx, ly);
    break;
  case 2:
    array_space::ij_to_r(neigh_ind, i, jp, lx, ly);
    break;
  case 3:
    array_space::ij_to_r(neigh_ind, im, j, lx, ly);
    break;
  case 4:
    array_space::ij_to_r(neigh_ind, im, jm, lx, ly);
    break;
  case 5:
    array_space::ij_to_r(neigh_ind, i, jm, lx, ly);
    break;
  }

  return neigh_ind;
}

int get_bond_triangular(const int site_1_ind, const int site_2_ind,
                        const int lx, const int ly) {

  for (int bond{0}; bond < 5; bond++) {
    int neighbour_of_bond{get_neighbour_triangular(site_1_ind, bond, lx, ly)};
    if (neighbour_of_bond == site_2_ind) {
      return bond;
    }
  }

  return 6;
}

const arr1i<6> &get_bond_permutation(int site_1_ind, int site_2_ind, int lx_m,
                               int ly_m) {
  int bond{get_bond_triangular(site_1_ind, site_2_ind, lx_m, ly_m)};
  return bond_permutation_struct::bond_permutations[static_cast<std::size_t>(
      bond)];
}

int get_face(const int orientation, const int bond) {
  std::size_t u_orientation{static_cast<std::size_t>(orientation)};
  const arr1i<6> &permutation_table{
      bond_permutation_struct::bond_permutations[u_orientation]};
  return permutation_table[static_cast<std::size_t>(bond)];
}

int get_interaction_index(const int site_1_orientation,
                          const int site_2_orientation, const int bond_1,
                          const int n_orientations) {
  int face_1{get_face(site_1_orientation, bond_1)};
  int bond_2{get_opposite_face(bond_1)};
  int face_2{get_face(site_2_orientation, bond_2)};

  return geometry_space::get_interaction_index(face_1, face_2, n_orientations);
}

int get_opposite_face(const int face) { return array_space::mod(face + 3, 6); }

} // namespace triangular_space
} // namespace geometry_space
