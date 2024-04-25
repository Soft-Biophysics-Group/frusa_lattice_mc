// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.
// Contents to be copy-pasted in app/main.cc to be ran

#include "model.h"
#include "mc_routines.h"
#include <iostream>


// Temporary includes for testing
#include "vector_utils.h"

int main(){

  geometry_space::Geometry geometry{};
  geometry_space::bond_struct triangular_bonds{geometry_space::triangular};
  std::cout << geometry << "\n\n" ;
  vec1i center_site_ijk {4, 0, 0};
  std::cout << "Creating site at coordinates " ;
  array_space::print_vector(std::cout, center_site_ijk);
  std::cout << "\n" ;
  int center_site_index {0};
  array_space::ijk_to_r(center_site_index, center_site_ijk[0],
      center_site_ijk[1], center_site_ijk[2], 5 , 5, 1);
  std::cout << "with index: " << center_site_index << '\n';
  int i_neigh {};
  int j_neigh {};
  int k_neigh {};
  int r_neigh {};
  for (int bond_index{0}; bond_index < 6; ++bond_index) {
    r_neigh = geometry.get_neighbour(center_site_index, bond_index);
    array_space::r_to_ijk(r_neigh, i_neigh, j_neigh, k_neigh, 5, 5, 1);
    std::cout << "Neighbour with bond " << bond_index << " has index "
              << r_neigh << " and coordinates {" << i_neigh << ", " << j_neigh
              << ", " << k_neigh << "}\n";
  }

  vec1i neigh_site_ijk  {1, 0, 0};
  int r_neigh_candidate {0};
  array_space::ijk_to_r(r_neigh_candidate, neigh_site_ijk[0], neigh_site_ijk[1],
                        neigh_site_ijk[2], 5, 5, 1);
  std::cout << r_neigh_candidate ;
  int bond_test {geometry.get_bond(center_site_index, r_neigh_candidate)};
  if (geometry.are_neighbours(bond_test)) {
    std::cout << "Center site with coordinates ";
    array_space::print_vector(std::cout, center_site_ijk);
    std::cout << " and neighbour site with coordinates " ;
    array_space::print_vector(std::cout, neigh_site_ijk);
    std::cout << " are neighbours with bond index " << bond_test << " corresponding to a vector " ;
    array_space::print_vector(
        std::cout,
        triangular_bonds.bond_array[static_cast<std::size_t>(bond_test)]);
    std::cout << '\n' ;
  } else {
    std::cout << "Center site with coordinates ";
    array_space::print_vector(std::cout, center_site_ijk);
    std::cout << " and neighbour site with coordinates " ;
    array_space::print_vector(std::cout, neigh_site_ijk);
    std::cout << "are not neigbours (bond index " << bond_test << ")\n" ;
  }

  int test_orientation {1};
  int test_bond{3};
  int test_coeff {geometry.get_interaction_coeff(test_orientation, test_bond)};
  std::cout << "Orientation " << test_orientation << " and bond " << test_bond
    << " are associated with matrix coefficient " << test_coeff << ".\n" ;

  return 0;
}
