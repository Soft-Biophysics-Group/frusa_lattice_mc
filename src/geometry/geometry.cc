#include "geometry.h"
#include "triangular.h"
#include "vector_utils.h"
#include <stdexcept>

using json = nlohmann::json;

namespace geometry_space {

bond_struct::bond_struct(lattice_options lattice) {
  switch (lattice) {
    case triangular:
      bond_permutation =
          triangular_space::bond_struct::bond_permutation;
      bond_array =
          triangular_space::bond_struct::bond_array;
      bond_index =
          triangular_space::bond_struct::bond_index;
      break;
    default:
      throw(std::runtime_error(
          "Wrong lattice option received duting bond_structure creation"));
  }
}

lattice_options get_lattice_from_str(std::string& lattice_str) {
  for (std::size_t i {0}; i < lattice_options::n_lattices; ++i) {
    if (lattice_str == lattice_str_arr[i]) {
      return static_cast<lattice_options>(i);
    }
  }
  throw(std::runtime_error("Invalid lattice name provided.\nOptions are: "
                           "chain, square, triangular, cubic, bcc, fcc"));
}

Geometry::Geometry(lattice_options lattice, int lx, int ly, int lz)
    : lattice_m{lattice}, lx_m{lx}, ly_m{ly}, lz_m{lz}, n_sites_m{lx * ly * lz},
      bond_struct_m{bond_struct(lattice)} {
  set_lattice_properties();
}

Geometry::Geometry(const std::string &geometry_input) {
  std::ifstream geometry_input_f{geometry_input};
  if (!geometry_input_f) {
    std::cerr << "Could not open geometry model parameters file" << '\n';
    exit(1);
  }

  json json_geometry = json::parse(geometry_input_f);

  std::string lattice_str{
      json_geometry["lattice_name"].template get<std::string>()};
  lattice_m = get_lattice_from_str(lattice_str);
  lx_m = json_geometry["lx"].template get<int>();
  ly_m = json_geometry["ly"].template get<int>();
  lz_m = json_geometry["lz"].template get<int>();
  n_sites_m = lx_m * ly_m * lz_m;
}

int Geometry::get_neighbour(const int site_ind, const int bond_ind) const {
  // Found at https://stackoverflow.com/a/26282004
    int i{};
    int j{};
    int k{};
    array_space::r_to_ijk(site_ind, i, j, k, lx_m, ly_m, lz_m);

    std::size_t u_bond_ind {static_cast<std::size_t>(bond_ind)};
    const vec1i& bond_direction {bond_struct_m.bond_array[u_bond_ind]};
    int i_neigh{i+bond_direction[0]};
    int j_neigh{j+bond_direction[1]};
    int k_neigh{k+bond_direction[2]};

    int neigh_ind{0};
    array_space::ijk_to_r(neigh_ind, i_neigh, j_neigh, k_neigh, lx_m, ly_m,
                          lz_m);
    return neigh_ind;
}

int Geometry::get_bond(const int site_1_ind, const int site_2_ind) const {
  // TODO Use the bond index array maybe?
  arr1i<3> site_1_ijk {0, 0, 0};
  array_space::r_to_ijk(site_1_ind, site_1_ijk[0], site_1_ijk[1], site_1_ijk[2],
                        lx_m, ly_m, lz_m);
  arr1i<3> site_2_to_site_1_ijk{0, 0, 0};
  array_space::r_to_ijk(site_2_ind, site_2_to_site_1_ijk[0],
                        site_2_to_site_1_ijk[1], site_2_to_site_1_ijk[2], lx_m,
                        ly_m, lz_m);
  for (std::size_t ind{0}; ind < 3; ind++) {
    site_2_to_site_1_ijk[ind] -= site_1_ijk[ind];
  }

  bool diff_found{bond_struct_m.bond_index.find(site_2_to_site_1_ijk) !=
                  bond_struct_m.bond_index.end()};
  if (diff_found) {
    return bond_struct_m.bond_index.at(site_2_to_site_1_ijk);
  }
  // Other implementation, based on comparing with neighbours
  // TODO See which one goes faster
  //for (int bond_index{0}; bond_index < n_neighbours_m; bond_index++) {
    //int neighbour_of_bond{get_neighbour(site_1_ind, bond_index)};
    //if (neighbour_of_bond == site_2_ind) {
      //return bond_index;
    //}
  //}
  // If not neighbours:
  return n_neighbours_m;
}

bool Geometry::are_neighbours(const int site_1_ind, const int site_2_ind) {
  // We have encoded two non-neghbour sites as having n_neighbours+1
  // orientations
  return (get_bond(site_1_ind, site_2_ind) != n_neighbours_m);
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

int Geometry::get_interaction_coeff(const int site_orientation, const int bond) const {
  std::size_t u_bond {static_cast<std::size_t>(bond)};
  const vec1i &permutation_table{bond_struct_m.bond_permutation[u_bond]};

  std::size_t u_orientation{static_cast<std::size_t>(site_orientation)};
  return permutation_table[u_orientation];
}

int Geometry::get_interaction_index(const int site_1_orientation,
                                    const int site_1_ind,
                                    const int site_2_orientation,
                                    const int site_2_ind) const {
  int bond{get_bond(site_1_ind, site_2_ind)};
  int coeff_1{get_interaction_coeff(site_1_orientation, bond)};
  int coeff_2{get_interaction_coeff(site_2_orientation, bond)};

  //Hash into a single coefficient for flattened vector
  int interaction_index{};
  array_space::ij_to_r(interaction_index, coeff_1, coeff_2, n_orientations_m, 1);
  return interaction_index;
}

double Geometry::get_interaction(const int site_1_orientation,
                                 const int site_1_ind,
                                 const int site_2_orientation,
                                 const int site_2_ind,
                                 const vec1d& flat_interaction_matrix) const {
  std::size_t interaction_index{static_cast<std::size_t>(get_interaction_index(
      site_1_orientation, site_1_ind, site_2_orientation, site_2_ind))};
  return flat_interaction_matrix[interaction_index];
}

std::ostream &operator<<(std::ostream &out, Geometry &geometry) {
  out << "Number of possible particle orientations:"
      << geometry.n_orientations_m << '\n';
  out << "Lattice dimensions:" << '(' << geometry.lx_m << ", " << geometry.ly_m
      << ", " << geometry.lz_m << ")\n";
  return out;
}

void Geometry::set_lattice_properties() {
  switch (lattice_m) {
  case lattice_options::triangular:
    n_neighbours_m = triangular_space::n_neighbours;
    n_orientations_m = triangular_space::n_orientations;
    break;
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}

} // namespace geometry_space
