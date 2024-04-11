#ifndef LATTICEPARTICLES_STATE_HEADER_H
#define LATTICEPARTICLES_STATE_HEADER_H

#include "lattice_particles_parameters.h"

#include <functional>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <sstream>

#include "vector_utils.h"

typedef std::mt19937 EngineType;
typedef std::uniform_int_distribution<int> int_dist;
typedef std::uniform_real_distribution<double> real_dist;

namespace lattice_particles_space{
  /*
   * Definitions required for the public routines of the model class
   */

  // Structure containing the characteristics of a lattice site.
  // Note that orientation 0 denotes an empty site, regardless of particle type
  class site_state {
    public:
      // TODO Do I need this constructor?
      site_state() = default;
      site_state(int type, int orientation, model_parameters_struct &params)
          : type_m{type}, orientation_m{orientation},
            n_orientations_m{params.n_orientations} {
          state_m = calc_state();
      }
      int get_type() { return type_m; };
      int get_orientation() { return orientation_m; } ;
      int get_state() { return state_m; };
      std::size_t get_site_index() { return site_index_m; }
      bool is_empty() { return orientation_m == 0; };
      void set_state(int type, int orientation) {
        type_m = type;
        orientation_m = orientation;
        state_m = calc_state();
      };
      void set_site_index(std::size_t new_index) { site_index_m = new_index; };
      friend std::ostream& operator<< (std::ostream& out, site_state &site);
      void swap_with(site_state& state2);
    private:
      int type_m {0};
      int orientation_m {0};
      int state_m {0};
      int n_orientations_m {0};
      std::size_t site_index_m {0};
      int calc_state() {
          return array_space::hash_into_state(type_m, orientation_m, n_orientations_m);
      };
  };

  using SiteVector = std::vector<site_state>;

  // Structure containing the characteristics of the state of the system:
  // n_types          - number of different particle trypes
  // n_orientations   - number of orientations a particle can take
  // n_states         - total number of states (type+orientation) a given
  //                    particle can take
  // lx, ly, lz       - dimensions of the lattice
  // n_sites          - total number of sites
  // n_particles      - number of particles of each type
  // lattice_sites    - array of each lattice state
  struct state_struct{
    int n_types {};
    int n_orientations {};
    int n_states {};
    int lx {};
    int ly {};
    int lz {};
    int n_sites {};
    vec1i n_particles {};
    SiteVector lattice_sites {};
    vec1i full_sites {};
    vec1i empty_sites {};
  };

  // Initialize the structural properties of the system, depending on the type
  // of parameters.initialize_option
  void initialize_state(state_struct &state,
                        model_parameters_struct &parameters);
  // Assign the indices to all the sites in lattice_sites
  void initialize_site_state_indices(SiteVector lattice_sites);
  // Print the current values of the structural properties of the system
  // TODO Implement this
  std::ostream& operator<< (std::ostream& out, state_struct &state);

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
  void initialize_state_from_file(state_struct &state,
                                  model_parameters_struct &parameters);

  // Initialize a state with a set random of particles uniformly distributed on
  // the lattice, with a given number of particles given in parameters
  void initialize_state_random_fixed_particle_numbers(
      state_struct &state, model_parameters_struct &parameters);
}

#endif
