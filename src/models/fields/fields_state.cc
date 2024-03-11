#include "fields_state.h" 

namespace fields_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_state(state_struct &state,
                        model_parameters_struct &parameters){
    
    state.ns = parameters.ns;
    state.Lx = parameters.Lx;
    state.Ly = parameters.Ly;
    state.Lz = parameters.Lz;
    state.Np = parameters.Np;
    
    // Set the total number of lattice sites and particle density
    state.N = state.Lx*state.Ly*state.Lz;
    state.rho_bar = (1.0*parameters.Np)/state.N;

    std::string option = parameters.initialize_option;

    try{
      if(option=="from_file"){
        initialize_state_from_file(state,parameters);
      }
      else if(option=="random"){
        initialize_state_random(state,parameters);
      }
      else if(option=="uniform"){
        initialize_state_uniform(state);
      }
      else{
        throw option;
      }
    }
    catch(std::string option){
      std::cout << "Incorrect initialization option: ''" << option << "''\n";
      exit(1);
    }

    for(int r=0;r<state.N;r++){
      if(state.local_density[r]>0){
        state.donor_list.push_back(r);
      }
      if(state.local_density[r]<1){
        state.acceptor_list.push_back(r);
      }
    }
  }

  void print_state(state_struct &state){
    std::cout << "\n---------------------------------------------\n";
    std::cout << "Current structural properties of the system\n";
    std::cout << "---------------------------------------------\n\n";
   
    std::cout << "Lattice dimensions:\n\n";
    std::cout << "Lx = " << state.Lx << "\n";
    std::cout << "Ly = " << state.Ly << "\n";
    std::cout << "Lz = " << state.Lz << "\n\n";

    std::cout << "Number of particles Np = " << state.Np << "\n\n";

    std::cout << "The total density of particles in the system is: ";
    std::cout << state.rho_bar << "\n\n";

    std::cout << "Current values of the fractional concentrations";
    std::cout << " and local densities:\n";
    for(int r=0;r<state.N;r++){

      int i,j,k;
      array_space::r_to_ijk(r,i,j,k,state.Lx,state.Ly,state.Lz);

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

  void save_state(state_struct &state, std::string state_output){
    
    std::ofstream state_f;
    state_f.open(state_output);
    if(!state_f){
      std::cerr << "Could not open "+state_output << std::endl;
      exit(1);
    }
    
    for(int r=0;r<state.N;r++){
      for(int s=0;s<state.ns;s++){
        state_f << std::setprecision(8) << state.concentration[r][s] << " ";
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
                                  model_parameters_struct &parameters){
    std::ifstream input_file;
    input_file.open(parameters.state_input);

    if(!input_file){
      std::cerr << "Unable to open file "+parameters.state_input;
      exit(1);
    }

    for(int r=0;r<state.N;r++){
      vec1d c_r;
      double rho_r = 0;
      for(int s=0;s<state.ns;s++){
        double c_r_s;
        input_file  >> c_r_s;
        c_r.push_back(c_r_s);
        rho_r+= c_r_s;
      }
      state.concentration.push_back(c_r);
      state.local_density.push_back(rho_r);
    }
  }

  void initialize_state_random(state_struct &state, 
                               model_parameters_struct &parameters){
    
    real_dist uniform_dist(0,1);

    for(int r=0;r<state.N;r++){
      
      vec1d c_r;
      double rho_r=0;
      
      for(int s=0;s<state.ns-1;s++){
        
        double c_r_s = uniform_dist(parameters.rng)*(state.rho_bar-rho_r);
        
        c_r.push_back(c_r_s);
        rho_r+= c_r_s;
      }

      c_r.push_back(state.rho_bar-rho_r);
      
      std::shuffle(std::begin(c_r), std::end(c_r), parameters.rng);

      state.concentration.push_back(c_r);
      state.local_density.push_back(state.rho_bar);
    }
  }

  void initialize_state_uniform(state_struct &state){
    
    for(int r=0;r<state.N;r++){
      
      vec1d c_r;
      
      for(int s=0;s<state.ns;s++){
        c_r.push_back(state.rho_bar/state.ns);
      }

      state.concentration.push_back(c_r);
      state.local_density.push_back(state.rho_bar);
    }
  }

  void update_state(int r, int index, int list_ind, double dc,
                    state_struct &state, double eps){

    state.concentration[r][index] += dc;
    state.local_density[r] += dc;

    if(std::find(state.donor_list.begin(),\
                 state.donor_list.end(),r)!=state.donor_list.end()){

      if(state.local_density[r]<eps){
        state.donor_list.erase(state.donor_list.begin()+list_ind);
      }
    }
    else{
      if(state.local_density[r]>eps){
        state.donor_list.push_back(r); 
      }
    }


    if(std::find(state.acceptor_list.begin(),\
                 state.acceptor_list.end(),r)!=state.acceptor_list.end()){

      if(state.local_density[r]>1-eps){
        state.acceptor_list.erase(state.acceptor_list.begin()+list_ind);
      }
    }
    else{
      if(state.local_density[r]<1-eps){
        state.acceptor_list.push_back(r); 
      }
    }
  }
}


