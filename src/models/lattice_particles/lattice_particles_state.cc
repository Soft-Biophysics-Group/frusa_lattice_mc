#include "lattice_particles_state.h"
#include <cstddef>
#include <ranges>

  /*
   * Definitions required for the public routines of the model class
   */
namespace lattice_particles_space{

  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_state(state_struct &state,
                        model_parameters_struct &parameters){

    state.n_types = parameters.n_types;
    state.n_orientations = parameters.n_orientations;
    state.lx = parameters.lx;
    state.ly = parameters.ly;
    state.lz = parameters.lz;
    state.n_sites = state.lx * state.ly * state.lz;
    state.n_particles = parameters.n_particles;
    state.lattice_sites =
        SiteVector(static_cast<std::size_t>(state.n_sites), site_state{});
    state.full_sites = vec1b(static_cast<std::size_t>(state.n_sites), false);

    // Set the total number of lattice sites and particle density
    state.n_sites = state.lx*state.ly*state.lz;

    std::string_view option = parameters.initialize_option;

    try{
      if(option=="from_file"){
        initialize_state_from_file(state,parameters);
      }
      else if(option=="random_fixed_particle_numbers"){
        // Segfault happens here
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

    // Start keeping track of full and empty sites
    for (std::size_t i {0}; i < state.lattice_sites.size(); i++){
      state.full_sites[i] = isempty(state.lattice_sites[i]);
    }
  }

  bool isempty(site_state& site) { return site.orientation == 0; };

  void save_state(state_struct &state, std::string state_output){
    std::ofstream state_f {state_output};
    if(!state_f){
      std::cerr << "Could not open "+state_output << std::endl;
      exit(1);
    }

    // Record the particle types on line 1
    for (std::size_t i = 0; i < state.lattice_sites.size(); i++) {
        state_f << state.lattice_sites[i].type << ' ';
    }
    state_f << '\n';
    // And the particle orientations on line 2
    for (std::size_t i = 0; i < state.lattice_sites.size(); i++) {
        state_f << state.lattice_sites[i].orientation << ' ';
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
                                  model_parameters_struct &parameters){
    std::ifstream input_file {parameters.state_input};

    if(!input_file){
      std::cerr << "Unable to open file "+parameters.state_input;
      exit(1);
    }

    // Fetching the orientations one by one
    std::string line {};
    int orientation {};
    std::getline(input_file, line);
    // TODO Would it work without the explicit definition? i.e. with:
    //while (std::stringstream(line) >> orientation) {...
    //Solution found at https://stackoverflow.com/a/20659156
    std::stringstream ss{line};
    while (ss >> orientation) {
        state.lattice_sites.push_back(site_state {0, orientation});
    }
    // Same for the particle types
    std::getline(input_file, line);
    int type {};
    ss = std::stringstream(line);
    // TODO This is a bit horrible. Find a better way.
    std::size_t site_index {0};
    while (ss >> type) {
      state.lattice_sites[site_index].type = type;
    }
  }

  void initialize_state_random_fixed_particle_numbers(
      state_struct &state, model_parameters_struct &parameters) {
      int_dist site_dist(0, state.n_sites);
      int_dist orientation_dist(0, parameters.n_orientations);

      // Fill the state until we get to the right number of particles of each
      // type
      for (std::size_t t {0}; t < static_cast<std::size_t>(parameters.n_types); t++)
      {
        for (std::size_t n {0}; n < static_cast<std::size_t>(parameters.n_particles[t]); n++) {
            std::size_t index{
                static_cast<std::size_t>(site_dist(parameters.rng))};
            int orientation{orientation_dist(parameters.rng)};
            state.lattice_sites[index] =
                site_state{orientation, static_cast<int>(t)};
        }
      }
  }

  //TODO I stopped here; continue
  void update_state(int index, int type, int orientation, state_struct &state) {
    // Unsigned index because the writes of this language hate you personnally
    // Yes, you
    std::size_t u_index {static_cast<std::size_t>(index)};
    state.lattice_sites[u_index].orientation = orientation;
    state.lattice_sites[u_index].type = type;
    // TODO LSP says ambiguous call below, but I don't see why.
    // See if compilation fails because of this.
    state.full_sites[u_index] = isempty(state.lattice_sites[u_index]);
  }

  void swap_sites(int index1, int index2, state_struct&state)  {
    std::size_t u_index1 {static_cast<std::size_t>(index1)};
    std::size_t u_index2 {static_cast<std::size_t>(index2)};
    site_state& site1 {state.lattice_sites[u_index1]};
    site_state& site2 {state.lattice_sites[u_index2]};

    int init_orientation_2 { site2.orientation };
    int init_type_2 { site2.type };

    // TODO Why is clangd throwing an ambiguous call error below?
    update_state(index2, site1.type, site1.orientation, state);
    update_state(index1, init_type_2, init_orientation_2, state);
  }
}
