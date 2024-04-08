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
            n_orientations{params.n_orientations} {
          state_m = calc_state();
      }
      int get_type() { return type_m; };
      int get_orientation() { return orientation_m; } ;
      int get_state() { return state_m; };
      bool is_empty() { return orientation_m == 0; };
      void set_state(int type, int orientation) {
        type_m = type;
        orientation_m = orientation;
        state_m = calc_state();
      };
      friend std::ostream& operator<< (std::ostream& out, site_state &site);
    private:
      int type_m {0};
      int orientation_m {0};
      int state_m {0};
      int n_orientations {0};
      int calc_state() {
          return array_space::hash_into_state(type_m, orientation_m, n_orientations);
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
  // full_sites       - Indices of booleans, where each one is true if the
  //                    corresponding site is true and empty otherwise
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
    vec1b full_sites {};
  };

  // Initialize the structural properties of the system, depending on the type
  // of parameters.initialize_option
  void initialize_state(state_struct &state,
                        model_parameters_struct &parameters);

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

  // Updates concentrations, local densities, and donor/acceptor lists
  // after a local update of a concentration field
  // index            - index of the lattice site to update
  // type             - new type of the particle at site index
  // orientation      - new orientation of the particle at site index
  // state            - state of the system before update
  void update_state(int index, int type, int orientation, state_struct &state);

  // Swap the states of sites at index1 and index2
  void swap_sites(int index1, int index2, state_struct& state);
}

#endif
