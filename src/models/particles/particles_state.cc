#include "particles_state.h"
#include "vector_utils.h"
#include <cstddef>
#include <numeric>
#include <sstream>

/*
 * Definitions required for the public routines of the model class
 */
namespace particles_space {

/*
 * Definitions required for the public routines of the model class
 */
void initialize_state(state_struct& state,
                      model_parameters_struct& parameters,
                      geometry_space::Geometry& geometry)
{
  state.n_types = parameters.n_types;
  state.n_orientations = geometry.get_n_orientations();
  state.n_particles = parameters.n_particles;
  state.n_states = state.n_types * state.n_orientations;
  state.n_sites = geometry.get_n_sites();
  std::string option = parameters.initialize_option;
  state.lattice_sites = SiteVector(option, state, parameters);
  state.full_empty_sites = FullEmptySites(state);
}

void save_state(state_struct& state, std::string& state_output)
{
  std::ofstream state_f;
  state_f.open(state_output);
  if (!state_f) {
    std::cerr << "Could not open " + state_output << std::endl;
    exit(1);
  }
  // Record the particle types on line 1
  for (int i {0}; i < state.n_sites; i++) {
    state_f << state.lattice_sites.get_type(i) << ' ';
  }

  state_f << '\n';
  // And the particle orientations on line 2
  for (int i {0}; i < state.n_sites; i++) {
    state_f << state.lattice_sites.get_orientation(i) << ' ';
  }
  state_f << '\n';
  state_f.close();
}
/*
 * End of the required definitions for the model class
 */

/*
 * Library-specific definitions
 */

/* -----------------------------
 * Functions to initialize state
 * -----------------------------*/

void initialize_state_from_file(vec1i& types,
                                vec1i& orientations,
                                model_parameters_struct& parameters)
{
  std::ifstream input_file {parameters.state_input};
  if (!input_file) {
    std::cerr << "Unable to open file " + parameters.state_input;
    exit(1);
  }
  // Fetching the orientations one by one
  std::string line {};
  int dummy {};
  std::getline(input_file, line);
  // Method found at https://stackoverflow.com/a/20659156 to convert text file
  // into int vector
  std::stringstream ss {line};
  while (ss >> dummy) {
    types.push_back(dummy);
  }
  // Same for the particle types
  std::getline(input_file, line);
  ss = std::stringstream(line);
  while (ss >> dummy) {
    orientations.push_back(dummy);
  }
  input_file.close();
}

void initialize_state_random_fixed_particle_numbers(
    vec1i& types,
    vec1i& orientations,
    state_struct& state,
    model_parameters_struct& parameters)
{
  int_dist orientation_dist(0, state.n_orientations - 1);
  // Fill the state until we get to the right number of particles of each
  // type
  std::size_t site_index {0};
  for (int type {0}; type < parameters.n_types; type++) {
    int n_particles_of_type {
        parameters.n_particles[static_cast<std::size_t>(type)]};
    for (int n {0}; n < n_particles_of_type; n++) {
      types.push_back(type);
      orientations.push_back(orientation_dist(parameters.rng));
      ++site_index;
    }
  }
  // Finish filling with empty sites
  while (site_index < static_cast<std::size_t>(state.n_sites)) {
    types.push_back(0);
    orientations.push_back(-1);
    site_index++;
  }
  // The lattice is now filled in order of particle types, so we need to
  // shuffle it
  // Since we have two arrays to shuffle, we first make a copy of our rng to
  // apply the same shuffle to each
  auto rng2 = parameters.rng;
  std::shuffle(types.begin(), types.end(), parameters.rng);
  std::shuffle(orientations.begin(), orientations.end(), rng2);
}

/* ----------------------------------
 *SiteVector class method definitions
 *----------------------------------*/

SiteVector::SiteVector(std::string& option,
                       state_struct& state,
                       model_parameters_struct& parameters)
    : n_orientations_m(state.n_orientations)
{
  try {
    if (option == "from_file") {
      initialize_state_from_file(types_m, orientations_m, parameters);
    } else if (option == "random") {
      initialize_state_random_fixed_particle_numbers(
          types_m, orientations_m, state, parameters);
    } else {
      throw option;
    }
  } catch (std::string option) {
    std::cout << "Incorrect initialization option: ''" << option << "''\n";
    exit(1);
  }
}

void print_state(state_struct& state)
{
  std::cout << state;
};

void SiteVector::swap_sites(const int index1, const int index2)
{
  std::size_t u_index1 {static_cast<std::size_t>(index1)};
  std::size_t u_index2 {static_cast<std::size_t>(index2)};

  std::swap(types_m[u_index1], types_m[u_index2]);
  std::swap(orientations_m[u_index1], orientations_m[u_index2]);
}

int SiteVector::get_state(const int site_index) const {
  // If site is empty, we'll return -1
  int state{-1};
  // Else we recycle Andrey's function to hash 2 integers into 1
  if (!is_empty(site_index)) {
    array_space::ij_to_r(state,
                         get_orientation(site_index),
                         get_type(site_index),
                         n_orientations_m,
                         1);
  }
  return state;
}

/* ---------------------------------------
 * FullEmptySites class method definitions
 * ---------------------------------------*/

FullEmptySites::FullEmptySites(state_struct& state)
{
  site_inds_to_full_empty_m = vec1s(static_cast<std::size_t>(state.n_sites));
  std::cout << state.n_sites << "\n";
  for (int i {0}; i < state.n_sites; i++) {
    std::size_t u_i {static_cast<std::size_t>(i)};
    if (state.lattice_sites.is_empty(i)) {
      empty_sites_indices_m.push_back(u_i);
      site_inds_to_full_empty_m[u_i] = empty_sites_indices_m.size() - 1;
    } else {
      full_sites_indices_m.push_back(u_i);
      site_inds_to_full_empty_m[u_i] = full_sites_indices_m.size() - 1;
    }
  }
}

int FullEmptySites::get_random_full_site(model_parameters_struct& parameters)
{
  int n_full_sites {static_cast<int>(full_sites_indices_m.size())};
  int_dist full_sites_dist(0, n_full_sites - 1);
  return static_cast<int>(full_sites_indices_m[static_cast<std::size_t>(
      full_sites_dist(parameters.rng))]);
}

int FullEmptySites::get_random_empty_site(model_parameters_struct& parameters)
{
  int n_empty_sites {static_cast<int>(empty_sites_indices_m.size())};
  int_dist empty_sites_dist(0, n_empty_sites - 1);
  return static_cast<int>(empty_sites_indices_m[static_cast<std::size_t>(
      empty_sites_dist(parameters.rng))]);
}

// TODO Check I got the right indices
void FullEmptySites::update_after_swap(const int initially_full_site,
                                       const int initially_empty_site)
{
  const std::size_t u_initially_empty_site {
      static_cast<std::size_t>(initially_empty_site)};
  const std::size_t u_initially_full_site {
      static_cast<std::size_t>(initially_full_site)};

  // Fetch the indices where empty and full arrays store full_index and
  // empty_index
  std::size_t index_in_empty_arr {
      site_inds_to_full_empty_m[u_initially_empty_site]};
  std::size_t index_in_full_arr {
      site_inds_to_full_empty_m[u_initially_full_site]};
  // Update the corresponding locations in the arrays of full and empty site
  // indices
  empty_sites_indices_m[index_in_empty_arr] = u_initially_full_site;
  full_sites_indices_m[index_in_full_arr] = u_initially_empty_site;
  site_inds_to_full_empty_m[u_initially_empty_site] = index_in_full_arr;
  site_inds_to_full_empty_m[u_initially_full_site] = index_in_empty_arr;
}

void swap_sites(state_struct& state, int site_1_index, int site_2_index)
{
  state.lattice_sites.swap_sites(site_1_index, site_2_index);
  if (state.lattice_sites.is_empty(site_1_index)) {
    // state.full_empty_sites.update_after_swap(site_2_index, site_1_index);
    state.full_empty_sites.update_after_swap(site_1_index, site_2_index);
  } else if (state.lattice_sites.is_empty(site_2_index)) {
    // state.full_empty_sites.update_after_swap(site_1_index, site_2_index);
    state.full_empty_sites.update_after_swap(site_2_index, site_1_index);
  }
}

std::ostream &operator<<(std::ostream &out, state_struct &state) {
  out << "Printing current system state\n";
  out << "Number of particle types: " << state.n_types << '\n';
  out << "Number of particles of each type: ";
  array_space::print_vector<int>(out, state.n_particles);
  out << '\n';
  out << "Particle type at each site:\n";
  array_space::print_vector(out, state.lattice_sites.types_m);
  out << '\n';
  out << "Particle orientation at each site:\n";
  array_space::print_vector(out, state.lattice_sites.orientations_m);
  out << '\n';
  out << "Indices of full sites:\n";
  array_space::print_vector(out, state.full_empty_sites.full_sites_indices_m);
  out << '\n';
  out << "Indices of empty sites:\n";
  array_space::print_vector(out, state.full_empty_sites.empty_sites_indices_m);
  out << '\n';
  out << "Mapping of lattice sites to full and empty lists:\n";
  array_space::print_vector(out,
                            state.full_empty_sites.site_inds_to_full_empty_m);
  out << '\n';
  return out;
}

} // namespace lattice_particles_space
