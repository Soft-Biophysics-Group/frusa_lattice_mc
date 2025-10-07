#ifndef GEOMETRY_H
#define GEOMETRY_H

/**
Andrey Zelenskiy, Vincent Ouazan-Reboul, 2024

This file contains class and functions needed to describe the geometry of
Bravais lattices on which we will be performing Monte-Carlo simulations, and in
particular the machinery for accesing the interactions between particles based
on which site they occupy and their relative orientations.

The list of lattices is contained in the lattice_options enum.

One key notion for the geometry part of the code is *bonds*.
A bond is a vector linking 2 neighbouring particles, given in lattice
coordinates and in the frame of reference of the lattice.
When we map an (orientation, orientation, bond) triplet to an interaction
coefficient, we first look at how each bond permutes the particle orientations,
or in other words: if I were to rotate the particle so that bond index b now
took the place of the bond 0, how would the orientation of the particle change?
This information is encoded in the bond_permutation array.
This is simple for 2D particles, where the faces are the only features that
permute under bond permutation. For 3D particles, this is trickier, and we
invite you to read the corresponding documentation.

Each lattice is associated with a `lattice.h` file defining its geometrical
properties. We invite you to visit these files for the lattices you will be
using in order to understand the conventions we use.
*/

#include "vector_utils.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <array>
#include <string>

using BondIndexMap = std::map<std::array<int,3> , int>;

namespace geometry_space
{
// Vectors of vectors giving the permutation indices for all the possible bond
// orientations.
enum lattice_options
{
  chain,
  square,
  triangular,
  cubic,
  bcc,
  fcc,
  n_lattices
};

static const inline std::array<std::string, lattice_options::n_lattices>
    lattice_str_arr {"chain", "square", "triangular", "cubic", "bcc", "fcc"};

lattice_options get_lattice_from_str(std::string& lattice_str);

/**
 * Structure containing the information about how neighbouring particles are
 * linked with each other
 */
struct bond_struct
{
  // TODO Write this constructor with switch cases
  bond_struct() = default;
  bond_struct(lattice_options lattice);
  // Which face we look at depending on the particle-particle bond and particle
  // orientation
  vec2i bond_permutation {};
  // Maps a bond index to the corresponding vector in lattice coordinates
  vec2i bond_array {};
  // And the reverse
  BondIndexMap bond_index {};
  // Maps a bond index to the opposite one.
  // Used for calculating particle contact indices: if particle 1 touches particle
  // 2 through bond, then particle 2 touches particle 1 through
  // opposite_bonds[bond]
  vec1i opposite_bonds {};
  // Number of neighbours per site
  int n_neighbours {};
};

std::ostream& operator<< (std::ostream& out, bond_struct& bonds);

/**
 * Class describing the geometry of a lattice and the particles occupying its
 * sites.
 * Typically built from an input file, which is by default
 * `input/model_params.json`.
 */
class Geometry {
public:
  // ----- CONSTRUCTORS -----
  Geometry(lattice_options lattice, int lx, int ly, int lz = 1);
  Geometry(const std::string& geometry_input = "input/model_params.json");

  // ----- SIMPLE GETTERS -----
  int get_n_orientations() const { return n_orientations_m; };
  int get_n_sites() const { return n_sites_m; };
  int get_n_neighbours() const { return n_neighbours_m; };

  // ----- GETTERS FOR NEIGHBOURING PARTICLES  AND SITES -----
  int get_neighbour(const int site_ind, const int bond_ind) const;
  int get_bond(const int site_1_ind, const int site_2_ind) const;
  int get_opposite_bond(const int bond) const
  {
    std::size_t u_bond {static_cast<std::size_t>(bond)};
    return bond_struct_m.opposite_bonds[u_bond];
  };
  bool are_neighbours(const int bond_index);
  bool are_neighbours(const int site_1_ind, const int site_2_ind);

  // ----- GETTERS FOR INTERACTIONS BETWEEN PARTICLES -----
  /**
   * Returns a one-particle index which will be hashed into a two-particles
   * contact index.
   * In 2D, this coefficient is simply the face that the particle is
   * presenting to the neighbour we are calculating the interaction with.
   * In 3D, it is a "generalized face" taking into account the overall
   * orientation of the particle.
   * Usually calculated as bond_type * (site_type + 1); runs from 1 to n_bonds
   * for particle type 1, bond_type * (site_type + 1) + 1 to 2 * (bond_type *
   * (site_type + 1)) for type 2, etc...
   * **/
  int get_interaction_coeff(const int site_orientation,
                            const int site_type,
                            const int bond) const;
  /**
   * Recovers the index of a particle-particle contact in the flattened
   * interaction matrix used as input.
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
  /**
   * Returns the value of the interacton energy between 2 particles sharing
   * neighbouring sites.
  */
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
  friend std::ostream& operator<<(std::ostream& out, Geometry& geometry);

private:
  // Which Bravais lattice we simulate
  lattice_options lattice_m {lattice_options::chain};
  // Lattice dimensions
  int lx_m {1};
  int ly_m {1};
  int lz_m {1};
  // Number of nerighbours per site
  int n_neighbours_m {2};
  // Number of orientations a particle can take
  int n_orientations_m {2};
  // Total number of sites on the lattice
  int n_sites_m {1};
  // Structure describing how neighbouring sites are linked
  bond_struct bond_struct_m {bond_struct(chain)};

  void set_lattice_properties();
};

}  // namespace geometry_space
#endif
