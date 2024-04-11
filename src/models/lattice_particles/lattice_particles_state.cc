#include "lattice_particles_state.h"
#include <cstddef>
#include <ranges>
#include "vector_utils.h"

/*
 * Definitions required for the public routines of the model class
 */
namespace lattice_particles_space{

  /*
   * Definitions required for the public routines of the model class
   */
void initialize_state(state_struct &state,
                      model_parameters_struct &parameters) {

  state.n_types = parameters.n_types;
  state.n_orientations = parameters.n_orientations;
  state.lx = parameters.lx;
  state.ly = parameters.ly;
  state.lz = parameters.lz;
  state.n_sites = state.lx * state.ly * state.lz;
  state.n_particles = parameters.n_particles;
  state.lattice_sites = SiteVector{};

  std::string option = parameters.initialize_option;

  try{
    if(option=="from_file"){
      initialize_state_from_file(state,parameters);
    }
    else if(option=="random_fixed_particle_numbers"){
      initialize_state_random_fixed_particle_numbers(state, parameters);
    }
    else{
      throw option;
    }
  }
  catch(std::string option){
    std::cout << "Incorrect initialization option: ''" << option << "''\n";
    exit(1);
  }
  initialize_site_state_indices(state.lattice_sites);
}

void site_state::swap_with(site_state& state2) {
  int temp_type {state2.type_m};
  int temp_orientation {state2.orientation_m};
  std::size_t temp_index {state2.site_index_m};

  state2.set_state(type_m, orientation_m);
  state2.set_site_index(site_index_m);
  set_state(temp_type, temp_orientation);
  set_site_index(state2.site_index_m);
}

void initialize_site_state_indices(SiteVector lattice_sites) {
   for (std::size_t i {0}; i<lattice_sites.size(); i++)
     lattice_sites[i].set_site_index(i);
}

void save_state(state_struct &state, std::string state_output) {

  std::ofstream state_f;
  state_f.open(state_output);
  if (!state_f) {
    std::cerr << "Could not open " + state_output << std::endl;
    exit(1);
  }

  // Record the particle types on line 1
  for (std::size_t i = 0; i < state.lattice_sites.size(); i++) {
    state_f << state.lattice_sites[i].get_type() << ' ';
  }
  state_f << '\n';
  // And the particle orientations on line 2
  for (std::size_t i = 0; i < state.lattice_sites.size(); i++) {
    state_f << state.lattice_sites[i].get_orientation() << ' ';
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
void initialize_state_from_file(state_struct &state,
                                model_parameters_struct &parameters) {
  std::ifstream input_file{parameters.state_input};

  if (!input_file) {
    std::cerr << "Unable to open file " + parameters.state_input;
    exit(1);
  }

  // Fetching the orientations one by one
  std::string line{};
  int type{};
  std::getline(input_file, line);
  // Solution found at https://stackoverflow.com/a/20659156
  std::stringstream ss{line};
  while (ss >> type) {
    state.lattice_sites.push_back(site_state{type, 0, parameters});
  }
  // Same for the particle types
  std::getline(input_file, line);
  int orientation{};
  std::stringstream ss2{line};
  // TODO This is a bit horrible. Find a better way.
  std::size_t site_index{0};
  while (ss2 >> orientation) {
    site_state &this_site{state.lattice_sites[site_index]};
    this_site.set_state(this_site.get_type(), orientation);
    ++site_index;
  }
  input_file.close();
}

std::ostream &operator<<(std::ostream &out, site_state &site) {
  out << "Index " << site.get_site_index() << " Type " << site.get_type() << ' '
      << " Orientation " << site.get_orientation() << " State " << site.get_state();
  return out;
}

void initialize_state_random_fixed_particle_numbers(
    state_struct &state, model_parameters_struct &parameters) {
  // Orientation has to start at 1 so that the site is not empty
  int_dist orientation_dist(1, parameters.n_orientations);
  // Fill the state until we get to the right number of particles of each
  // type
  std::size_t site_index{0};
  for (int type{0}; type < parameters.n_types; type++) {
    int n_particles_of_type{
        parameters.n_particles[static_cast<std::size_t>(type)]};
    for (int n{0}; n < n_particles_of_type; n++) {
      state.lattice_sites.push_back(
          site_state(type, orientation_dist(parameters.rng), parameters));
      ++site_index;
    }
  }
  // Finish filling with empty sites
  while (site_index < static_cast<std::size_t>(state.n_sites)) {
    state.lattice_sites.push_back(site_state(0, 0, parameters));
    site_index++;
  }
  // The lattice is now filled in order of particle types, so we need to
  // shuffle it
  std::shuffle(state.lattice_sites.begin(), state.lattice_sites.end(),
               parameters.rng);
}

std::ostream& operator<< (std::ostream& out, state_struct &state) {
  out << "Printing current system state\n" ;
  out << "Number of particle types: " << state.n_types << '\n' ;
  out << "Number of particle orientations: " << state.n_orientations << '\n' ;
  out << "Lattice dimensions: (" << state.lx << ", " << state.ly << ", "
      << state.lz << ")\n";
  out << "Number of lattice sites: " << state.n_sites << '\n' ;
  out << "Number of particles of each type: ";
  array_space::print_vector<int>(out, state.n_particles);
  out << '\n';
  out << "Printing site states:";
  for (auto& site : state.lattice_sites)
    out << site;

  return out;
}

} // namespace lattice_particles_space
