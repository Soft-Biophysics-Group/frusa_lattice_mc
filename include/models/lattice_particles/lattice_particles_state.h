#ifndef LATTICEPARTICLES_STATE_HEADER_H
#define LATTICEPARTICLES_STATE_HEADER_H

#include "lattice_particles_parameters.h"

#include <functional>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>

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
  struct site_state {
    int type {0};
    int orientation {0};
  };
  bool isempty(site_state& site);
  std::ostream& operator<< (std::ostream& out, site_state &site);

  using SiteVector = std::vector<site_state>;
  using SiteBool = std::pair<site_state, bool>;

  // Structure containing the characteristics of the state of the system:
  // n_types          - number of different particle trypes
  // n_orientations   - number of orientations a particle can take
  // lx, ly, lz       - dimensions of the lattice
  // n_sites          - total number of sites
  // n_particles      - number of particles
  // lattice_sites    - array of each lattice state
  // full_sites       - Indices of all sites containing a particle
  // empty_sites      - Indices of all empty sites
  struct state_struct{
    int n_types {};
    int n_orientations {};
    int lx {};
    int ly {};
    int lz {};
    int n_sites {};
    int n_particles {};
    SiteVector lattice_sites {};
    vec1i full_sites {};
    vec1i empty_sites {};
  };

  // Initialize the structural properties of the system, depending on the type
  // of parameters.initialize_option
  void initialize_state(state_struct &state,
                        model_parameters_struct &parameters);

  // Print the current values of the structural properties of the system
  void print_state(state_struct &state);

  // Save the fractional concentrations to a file "state_output"
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

  void initialize_state_random(state_struct &state,
                               model_parameters_struct &parameters);

  void initialize_state_uniform(state_struct &state);

  // Updates concentrations, local densities, and donor/acceptor lists
  // after a local update of a concentration field
  // r        - lattice position of the updated field
  // index    - component of the updated field
  // list_ind - index of the updated field in either the donor or acceptor list
  // dc       - amount by which the value of the field is changed
  // state    - state of the system before update
  void update_state(int r, int index, int list_ind, double dc,
                    state_struct &state, double eps=1e-8);
}

#endif
