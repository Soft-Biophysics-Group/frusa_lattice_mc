#include "fields_state.h" 

namespace fields_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_state(state_struct &state,
                        model_parameters_struct parameters){
    
    state.ns = parameters.ns;
    state.Lx = parameters.Lx;
    state.Ly = parameters.Ly;
    state.Lz = parameters.Lz;
    
    // Set the total number of lattice sites and particle density
    state.N = state.Lx*state.Ly*state.Lz;
    state.rho_bar = parameters.Np/state.N;

    std::string option = parameters.initialize_option;

    try{
      if(option=="from_file"){
        initialize_state_from_file(state,parameters);
      }
      else if(option=="random"){
        initialize_state_random(state,parameters);
      }
      else if(option=="uniform"){
        initialize_state_uniform(state,parameters);
      }
      else{
        throw option;
      }
    }
    catch(std::string option){
      std::cout << "Incorrect initialization option: ''" << option << "''\n";
      exit(1);
    }
  }

  void print_state(state_struct state, model_parameters_struct parameters){
    std::cout << "\n---------------------------------------------\n";
    std::cout << "Current structural properties of the system\n";
    std::cout << "---------------------------------------------\n\n";
    
    std::cout << "The total density of particles in the system is:";
    std::cout << state.rho_bar << "\n\n";

    std::cout << "Current values of the fractional concentrations";
    std::cout << " and local densities:\n";
    for(int r=0;r<state.N;r++){

      //TODO define get_indices in utils
      int i,j,k;
      r_to_indices(r,i,j,k);

      std::cout << "i = " << i << " ";
      std::cout << "j = " << j << " ";
      std::cout << "k = " << k << " ";
      
      std::cout << "c(r) = " << " ";

      for(int s=0;s<state.ns;s++){
        std::cout << state.concentration[r][s] << " ";
      }

      std::cout << " rho(r) = " << " " << state.local_density[r];
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }

  void save_state(state_struct state, model_parameters_struct parameters){
    
    std::string output_address = parameters.output_address;

    std::ofstream state_f;
    state_f.open(output_address);
    if(!state_f){
      std::cerr << "Could not open "+output_address << std::endl;
      exit(1);
    }
    
    for(int r=0;r<state.N;r++){
      for(int s=0;s<state.ns;s++){
        state_f << std::setprecision(8) << concentrations[r][s] << " ";
      }
      state_f << "\n";
    }
    state_f.close();
  }
  /* 
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
  void initialize_state_from_file(state_struct &state, 
                                  model_parameters_struct parameters){

  }

  void initialize_state_random(state_struct &state, 
                               model_parameters_struct parameters){
    
    uniform_dist  = real_dist(0,1);
    rng = parameters.rng;

    for(int r=0;r<state.N;r++){
      
      vec1d c_r;
      double rho_r;
      
      for(int s=0;s<state.ns-1;s++){
        
        double c_r_s = uniform_dist(rng)*rho_bar;
        
        c_r.push_back(c_r_s);
        rho_r+= c_r_s;
      }

      c_r.push_back(state.rho_bar-rho_r);

      state.concentrations.push_back(c_r);
      state.local_density.push_back(state.rho_bar);
    }
  }

  void initialize_state_uniform(state_struct &state){
    
    for(int r=0;r<state.N;r++){
      
      vec1d c_r;
      
      for(int s=0;s<state.ns;s++){
        c_r.push_back(state.rho_bar/state.ns);
      }

      state.concentrations.push_back(c_r);
      state.local_density.push_back(state.rho_bar);
    }
  }

  vec1i get_neighbours(int r, state_struct state){
    //TODO redefine for general Lx Ly Lz

    vec1i neighbours;

    int rm = (N+((r-1)%N))%N;
    int rp = (r+1)%N;

    neighbours.push_back({0,rm});
    neighbours.push_back({1,rp});

    return neighbours;
  }
}


