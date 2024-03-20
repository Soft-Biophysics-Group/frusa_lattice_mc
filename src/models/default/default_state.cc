#include "default_state.h" 

namespace default_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_state(state_struct &state,
                        model_parameters_struct &parameters){
    
    state.N = parameters.N;
    
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
  }

  void print_state(state_struct &state){
    std::cout << "\n---------------------------------------------\n";
    std::cout << "Current state of the system\n";
    std::cout << "---------------------------------------------\n\n";
   
    std::cout << "Number of particles N = " << state.N << "\n\n";

    std::cout << "The average occupation number is: ";
    std::cout << state.average_occupation << "\n\n";

    std::cout << "Current state occupation:\n";
    for(int r=0;r<state.N;r++){

      std::cout << "r = " << r << " ";
      std::cout << "n(r) = " << " ";

      std::cout << state.occupation[r] << "\n";
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
      state_f << state.occupation[r] << "\n";
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

    double n_av = 0;
    for(int r=0;r<state.N;r++){
      double n_r = 0;
      input_file >> n_r;
      
      state.occupation.push_back(n_r);

      n_av+= n_r;  
    }
    
    state.average_occupation = n_av/state.N;
  }

  void initialize_state_random(state_struct &state, 
                               model_parameters_struct &parameters){
    
    int_dist binary_dist(0,1);

    double n_av = 0;

    for(int r=0;r<state.N;r++){
      
      int n_r = binary_dist(parameters.rng);
        
      state.occupation.push_back(n_r);

      n_av+= n_r;  
    }
    
    state.average_occupation = n_av/state.N;
  }

  void initialize_state_uniform(state_struct &state){
    
    for(int r=0;r<state.N;r++){
      state.occupation.push_back(1);
    }
    
    state.average_occupation = 1;
  }
}


