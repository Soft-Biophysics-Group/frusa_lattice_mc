#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "vector_utils.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <array>
#include <string>

using BondIndexMap = std::map<std::array<int,3> , int>;

namespace geometry_space {
// Vectors of vectors giving the permutation indices for all the possible bond
// orientations.
enum lattice_options {
  chain,
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

struct bond_struct{
  // TODO Write this constructor with switch cases
  bond_struct() = default;
  bond_struct(lattice_options lattice);
  vec2i bond_permutation{};
  vec2i bond_array{};
  BondIndexMap bond_index{};
};

std::ostream& operator<< (std::ostream& out, bond_struct& bonds);

class Geometry {
public:
  Geometry(lattice_options lattice, int lx, int ly, int lz = 1);
  Geometry(const std::string &geometry_input = "input/model_params.json");
  // Simple getters
  int get_n_orientations() const { return n_orientations_m; };
  int get_n_sites() const { return n_sites_m; };
  int get_n_neighbours() const { return n_neighbours_m; };
  // More involved functions
  int get_neighbour(const int site_ind, const int bond_ind) const;
  int get_bond(const int site_1_ind, const int site_2_ind) const;
  // template <int N>
  // arr1i<N>& get_bond_permutation(const int site_1_ind, const int site_2_ind)
  // const;
  int get_interaction_coeff(const int site_orientation, const int site_type,
                            const int bond) const;
  /**
   * Return the index of a contact between two particles in the flattened
   * interactions array.
  * **/
  int get_interaction_index(const int site_1_orientation,
                            const int site_1_type,
                            const int site_2_orientation,
                            const int site_2_type,
                            const int bond,
                            const int n_types) const;
  int get_interaction_index(const int site_1_orientation,
                            const int site_1_type,
                            const int site_1_ind,
                            const int site_2_orientation,
                            const int site_2_type,
                            const int site_2_ind,
                            const int n_types) const
  {
    return get_interaction_index(site_1_orientation,
                                 site_1_type,
                                 site_2_orientation,
                                 site_2_type,
                                 get_bond(site_1_ind, site_2_ind),
                                 n_types);
  };
  double get_interaction(const int site_1_orientation,
                         const int site_1_type,
                         const int site_2_orientation,
                         const int site_2_type,
                         const int bond,
                         const int n_types,
                         const vec1d& flat_interaction_matrix) const;
  double get_interaction(const int site_1_orientation,
                         const int site_1_type,
                         const int site_1_ind,
                         const int site_2_orientation,
                         const int site_2_type,
                         const int site_2_ind,
                         const int n_types,
                         const vec1d& flat_interaction_matrix) const
  {
    return get_interaction(site_1_orientation,
                           site_1_type,
                           site_2_orientation,
                           site_2_type,
                           get_bond(site_1_ind, site_2_ind),
                           n_types,
                           flat_interaction_matrix);
  }
  bool are_neighbours(const int bond_index);
  bool are_neighbours(const int site_1_ind, const int site_2_ind);
  friend std::ostream& operator<<(std::ostream& out, Geometry& geometry);

private:
  lattice_options lattice_m {lattice_options::chain};
  int lx_m {1};
  int ly_m {1};
  int lz_m {1};
  int n_neighbours_m {2};
  int n_orientations_m {2};
  int n_sites_m {1};
  bond_struct bond_struct_m;
  void set_lattice_properties();
};

} // namespace geometry_space
#endif
