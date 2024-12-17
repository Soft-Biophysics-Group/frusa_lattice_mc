#ifndef PARTICLES_STATE_HEADER_H
#define PARTICLES_STATE_HEADER_H

#include "particles_parameters.h"
#include "geometry.h"

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

namespace particles_space
{

// Forward declaration for member classes
struct state_struct;

/*
 * Definitions required for the public routines of the model class
 */

/*
 * Class that stores the state of the full lattice, as two lists: one of
 * particle types and the other of particle orientations. Orientation -1 will
 * always be empty. Created along with the state according to one of the options
 * which must be supplied in the config file.
 */
class SiteVector {
public:
  SiteVector() = default;
  /**
   * Parameters:
   * - option: string which can take 2 values, "from_file", in which case state
   *   is initialized from the contents of a file with appropriate format, and
   *   "random", in which case the state is randomly initialized.
   * - state: state_struct class instance
   * - parameters: appropriate model_parameters_struct instance
   * **/
  SiteVector(std::string& option,
             state_struct& state,
             model_parameters_struct& parameters);

  // Simple getters and setters
  int get_type(const int site_index) const
  {
    return types_m[static_cast<std::size_t>(site_index)];
  };
  int get_orientation(const int site_index) const
  {
    return orientations_m[static_cast<std::size_t>(site_index)];
  };
  bool is_empty(const int site_index) const
  {
    return orientations_m[static_cast<std::size_t>(site_index)] == -1;
  };
  int get_state(const int site_index) const;
  void set_type(const int site_index, const int new_type)
  {
    types_m[static_cast<std::size_t>(site_index)] = new_type;
  }
  void set_orientation(const int site_index, const int new_orientation)
  {
    orientations_m[static_cast<std::size_t>(site_index)] = new_orientation;
  }
  void set_site(const int site_index,
                const int new_type,
                const int new_orientation)
  {
    set_type(site_index, new_type);
    set_orientation(site_index, new_orientation);
  }

  /**
   * Swap contents (particle orientations and types) of sites with indices
   * `index1` and `index2`
   * **/
  void swap_sites(const int index1, const int index2);
  friend std::ostream& operator<<(std::ostream& out, state_struct& state);

private:
  // Vector of particle types
  vec1i types_m{};
  // Vector of particle orientations. -1 is empty
  vec1i orientations_m{};
  // Number of orientations our particles can take
  int n_orientations_m{};
};

/*
 * Class which keeps track of which lattice sites are full and empty.
 * To be updated whenever we swap full and empty sites.
 * Useful to pick a full/empty site at random.
 */
class FullEmptySites
{
public:
  FullEmptySites() = default;
  FullEmptySites(state_struct &state);
  void update_after_swap(const int initially_full_index,
                         const int initially_empty_index);
  // Simple getters
  int get_n_full_sites()
  {
    return static_cast<int>(full_sites_indices_m.size());
  };
  int get_n_empty_sites()
  {
    return static_cast<int>(empty_sites_indices_m.size());
  };
  // Get random sites for MC moves
  int get_random_full_site(model_parameters_struct& parameters);
  int get_random_empty_site(model_parameters_struct& parameters);
  friend std::ostream& operator<<(std::ostream& out, state_struct& state);

private:
  // Vector of full site indices
  vec1s full_sites_indices_m{};
  // Vector of empty site indices
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
// lattice_sites    - See SiteVector class
// full_empty_sites - See FullEmptySites class
struct state_struct {
  // Number of particle types
  int n_types{};
  // Number of orietations a particle can take
  int n_orientations{};
  // Total number of states a particle can take
  int n_states{};
  // Number of lattice sites
  int n_sites{};
  // Vector of the number of particles of each type
  vec1i n_particles{};
  // Class describing the state of each lattice site
  SiteVector lattice_sites{};
  // Class keeping track of which sites are full and which sites are empty
  FullEmptySites full_empty_sites{};
};

// Initialize the structural properties of the system, depending on the type
// of parameters.initialize_option
void initialize_state(state_struct& state,
                      model_parameters_struct& parameters,
                      geometry_space::Geometry& geometry);

// Print the current values of the structural properties of the system
std::ostream& operator<<(std::ostream& out, state_struct& state);
// Thin wrapper around operator<< to make this header consistent with the rest
void print_state(state_struct &state);

// Save the state of the lattice to a file "state_output"
// First line of the file is a list of all the sites' particle types
// Second line is the list of particle orientations
void save_state(state_struct& state, std::string& state_output);

/*
 * End of the required definitions for the model class
 */

/*
 * Library-specific definitions
 */

// Various methods for state initialization
// Initialize state from a file specified in the model_parameters_struct
// object. Will fill up the types and orientations array of a SiteVector class.
void initialize_state_from_file(vec1i& types,
                                vec1i& orientations,
                                model_parameters_struct& parameters);

// Initialize a state with a set random of particles uniformly distributed on
// the lattice, with a given number of particles given in parameters
// types and orientations are the arrays being filled
void initialize_state_random_fixed_particle_numbers(
    vec1i& types,
    vec1i& orientations,
    state_struct& state,
    model_parameters_struct& parameters);

// Exchange the states of site_1 and site_2, updating both the SiteVector and
// FullEmptySites objects
void swap_sites(state_struct& state, int site_1_index, int site_2_index);

}  // namespace particles_space
#endif
