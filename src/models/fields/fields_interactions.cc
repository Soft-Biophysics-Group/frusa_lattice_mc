#include "fields_interactions.h"

namespace fields_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void initialize_interactions(state_struct &state, 
                               interactions_struct &interactions,
                               model_parameters_struct &parameters){

    // Initialize model parameters from user input file
    interactions.coupling_matrix = get_coupling_matrix(parameters.couplings);
    interactions.T_model         = parameters.couplings[0];

    // Calculate mean-field free energy of the initial state
    interactions.energy      = get_energy(state, interactions.coupling_matrix);
    interactions.entropy     = get_entropy(state, interactions.T_model);
    interactions.free_energy = interactions.energy+interactions.entropy;
  }

  void print_interactions(state_struct &state,
                          interactions_struct &interactions){
    std::cout << "\n--------------------------------------------\n";
    std::cout << "Current interaction properties of the system\n";
    std::cout << "--------------------------------------------\n\n";
    
    std::cout << "The coupling matrix on the bond along the x-axis is:\n\n";
    for(int i=0;i<interactions.coupling_matrix.size();i++){
      for(int j=0;j<interactions.coupling_matrix[i].size();j++){
        std::cout << interactions.coupling_matrix[0][i][j] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Mean-field temperature T_model = ";
    std::cout << interactions.T_model << "\n\n";

    std::cout << "Current value mean-field energy, entropy, and";
    std::cout << " free energy per particle:\n\n";
    std::cout << "energy = " << interactions.energy/state.Np << "\n";
    std::cout << "entropy = " << interactions.entropy/state.Np << "\n";
    std::cout << "free energy = " << interactions.free_energy/state.Np << "\n";
    std::cout << "\n";
  }
  /* 
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
  double get_energy(state_struct &state, vec3d coupling_matrix){

    double energy = 0;

    for(int r1=0;r1<state.N;r1++){

      vec1d c_r1 = state.concentration[r1];
      
      vec1i n = get_neighbours(r1,state.Lx,state.Ly,state.Lz);

      for(int k=0;k<n.size();k++){
        int r2 = n[k];
        vec1d c_r2= state.concentration[r2];

        for(int a=0;a<state.ns;a++){
          for(int b=0;b<state.ns;b++){
            energy += coupling_matrix[k][a][b]*c_r1[a]*c_r2[b];
          }
        }
      }
    }
    return energy/2;
  }

  double get_entropy(state_struct &state, double T_model, double eps){
    
    double entropy = 0;

    for(int r=0;r<state.N;r++){
      if(1-state.local_density[r]>eps){
        entropy += T_model*(1-state.local_density[r])*\
                   log(1-state.local_density[r]);
      }
      for(int a=0;a<state.ns;a++){
        if(state.concentration[r][a]>eps){
          entropy += T_model*state.concentration[r][a]*\
                     log(state.concentration[r][a]);
        }
      }
    }
    return entropy;
  }
}
