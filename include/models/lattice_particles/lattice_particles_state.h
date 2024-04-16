#ifndef LATTICEPARTICLES_STATE_HEADER_H
#define LATTICEPARTICLES_STATE_HEADER_H

#include "lattice_particles_parameters.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "vector_utils.h"

typedef std::mt19937 EngineType;
typedef std::uniform_int_distribution<int> int_dist;
typedef std::uniform_real_distribution<double> real_dist;

namespace lattice_particles_space {
// Forward declaration for member classes
struct state_struct;

/*
 * Definitions required for the public routines of the model class
 */

// Get a hashed state integer from particle characteristics, returning 0 if
// the site is empty
int hash_into_state(int type, int orientation, int n_orientations);

/*
 * Class that stores the state of the full lattice.
 * Created along with the state according to one of the options which must be
 * supplied in the config file.
 */
class SiteVector {
public:
  SiteVector() = default;
  SiteVector(std::string &option, state_struct &state,
             model_parameters_struct &parameters);
  int get_type(const int site_index) const {
    return types_m[static_cast<std::size_t>(site_index)];
  };
  int get_orientation(const int site_index) const {
    return orientations_m[static_cast<std::size_t>(site_index)];
  };
  bool is_empty(const int site_index) const {
    return orientations_m[static_cast<std::size_t>(site_index)] == 0;
  };
  int get_state(const int site_index) const {
    return hash_into_state(
        types_m[static_cast<std::size_t>(site_index)],
        orientations_m[static_cast<std::size_t>(site_index)], n_orientations_m);
  }
  void set_type(const int site_index, const int new_type) {
    types_m[static_cast<std::size_t>(site_index)] = new_type;
  }
  void set_orientation(const int site_index, const int new_orientation) {
    orientations_m[static_cast<std::size_t>(site_index)] = new_orientation;
  }
  void set_site(const int site_index, const int new_type,
                const int new_orientation) {
    set_type(site_index, new_type);
    set_orientation(site_index, new_orientation);
  }
  void swap_sites(const int index1, const int index2);
  friend std::ostream &operator<<(std::ostream &out, state_struct &state);

private:
  vec1i types_m{};
  vec1i orientations_m{};
  int n_orientations_m{};
};

/*
 * Class which keeps track of which lattice sites are full and empty.
 * To be updated whenever we swap full and empty sites.
 */
class FullEmptySites {
public:
  FullEmptySites() = default;
  FullEmptySites(state_struct &state);
  void update_after_swap(const int initially_full_index,
                         const int initially_empty_index);
  int get_n_full_sites() {
    return static_cast<int>(full_sites_indices_m.size());
  };
  int get_n_empty_sites() {
    return static_cast<int>(empty_sites_indices_m.size());
  };
  int get_random_full_site(model_parameters_struct& parameters);
  int get_random_empty_site(model_parameters_struct& parameters);
  friend std::ostream &operator<<(std::ostream &out, state_struct &state);

private:
  vec1s full_sites_indices_m{};
  vec1s empty_sites_indices_m{};
  // map of each site index to the corresponding coefficient in the
  // full_/empty_sites arrays
  vec1s site_inds_to_full_empty_m{};
};

// Structure containing the characteristics of the state of the system:
// n_types          - number of different particle trypes
// n_orientations   - number of orientations a particle can take
// n_states         - total number of states (type+orientation) a given
//                    particle can take
// lx, ly, lz       - dimensions of the lattice
// n_sites          - total number of sites
// n_particles      - number of particles of each type
// lattice_sites    - array of each lattice state
struct state_struct {
  int n_types{};
  int n_orientations{};
  int n_states{};
  int lx{};
  int ly{};
  int lz{};
  int n_sites{};
  vec1i n_particles{};
  SiteVector lattice_sites{};
  FullEmptySites full_empty_sites{};
};

// Initialize the structural properties of the system, depending on the type
// of parameters.initialize_option
void initialize_state(state_struct &state, model_parameters_struct &parameters);
void fill_full_empty_arrays(state_struct &state);

// Print the current values of the structural properties of the system
// TODO Implement this
std::ostream &operator<<(std::ostream &out, state_struct &state);

// Save the fractional concentrations to a file "state_output"
// First line of the file is a list of all the sites' particle types
// Second line is the list of particle orientations
void save_state(state_struct &state, std::string state_output);

/*
 * End of the required definitions for the model class
 */

/*
 * Library-specific definitions
 */

// Various methods for state initialization
void initialize_state_from_file(vec1i& types, vec1i& orientations,
                                model_parameters_struct &parameters);

// Initialize a state with a set random of particles uniformly distributed on
// the lattice, with a given number of particles given in parameters
void initialize_state_random_fixed_particle_numbers(
    vec1i& types, vec1i& orientations, state_struct &state,
    model_parameters_struct &parameters);

void swap_sites(state_struct &state, int site_1_index, int site_2_index);

int type_of_state(int state, int n_orientations);

void print_state(state_struct& state);

} // namespace lattice_particles_space
#endif
