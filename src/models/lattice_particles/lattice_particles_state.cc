#include "lattice_particles_state.h"

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
    state.n_particles = parameters.n_particles;

    // Set the total number of lattice sites and particle density
    state.n_sites = state.lx*state.ly*state.lz;

    std::string_view option = parameters.initialize_option;

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

    // Start keeping track of full and empty sites
    for (int i {0}; i < state.n_sites; i++){
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
    for (int i = 0; i < state.n_sites; i++) {
        state_f << state.lattice_sites[i].type << ' ';
    }
    state_f << '\n';
    // And the particle orientations on line 2
    for (int i = 0; i < state.n_sites; i++) {
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

    SiteVector lattice_sites {};
    // Fetching the orientations one by one
    std::string line {};
    int orientation {};
    std::getline(input_file, line);
    // TODO Would it work without the explicit definition? i.e. with:
    //while (std::stringstream(line) >> orientation) {...
    //Solution found at https://stackoverflow.com/a/20659156
    std::stringstream ss{line};
    while (ss >> orientation) {
        lattice_sites.push_back(site_state {0, orientation});
    }
    // Same for the particle types
    std::getline(input_file, line);
    int type {};
    ss = std::stringstream(line);
    // TODO This is a bit horrible. Find a better way.
    int site_index {0};
    while (ss >> type) {
      lattice_sites[site_index].type = type;
    }
  }

  void initialize_state_random_fixed_particle_numbers(
      state_struct &state, model_parameters_struct &parameters) {
      int_dist site_dist(0, state.n_sites);
      int_dist orientation_dist(0, parameters.n_orientations);

      // Fill the state until we get to the right number of particles of each
      // type
      for (int t {0}; t < parameters.n_types; t++)
      {
        for (int n {0}; n < parameters.n_particles[n]; n++) {
          int index{site_dist(parameters.rng)};
          int orientation{orientation_dist(parameters.rng)};
          state.lattice_sites[index] = site_state{orientation, t};
        }
      }
  }

  //TODO I stopped here; continue
  void update_state(int index, int type, int orientation, state_struct &state) {
    state.lattice_sites[index].orientation = orientation;
    state.lattice_sites[index].type = type;
    // TODO LSP says ambiguous call below, but I don't see why.
    // See if compilation fails because of this.
    state.full_sites[index] = isempty(state.lattice_sites[index]);
  }

  void swap_sites(int index1, int index2, state_struct&state)  {
    site_state& site1 {state.lattice_sites[index1]};
    site_state& site2 {state.lattice_sites[index2]};

    int init_orientation_2 { site2.orientation };
    int init_type_2 { site2.type };

    // TODO Why is clangd throwing an ambiguous call error below?
    update_state(index2, site1.type, site1.orientation, state);
    update_state(index1, init_type_2, init_orientation_2, state);
  }
}
