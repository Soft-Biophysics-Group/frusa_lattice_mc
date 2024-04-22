#include "geometry.h"
#include "triangular.h"
#include <stdexcept>

using json = nlohmann::json;

namespace geometry_space {

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
    : lattice_m{lattice}, lx_m{lx}, ly_m{ly}, lz_m{lz} {
  set_lattice_properties();
}

Geometry::Geometry(std::string& geometry_input)
{
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

void Geometry::set_lattice_properties() {
  switch (lattice_m) {
  case lattice_options::triangular:
    n_neighbours_m = 6;
    n_orientations_m = 6;
    break;
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
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
