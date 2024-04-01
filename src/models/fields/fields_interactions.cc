// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

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

  void print_energy(state_struct &state, interactions_struct &interactions){
    std::cout << "energy = " << interactions.energy/state.Np << "\n";
    std::cout << "entropy = " << interactions.entropy/state.Np << "\n";
    std::cout << "free energy = " << interactions.free_energy/state.Np << "\n\n";
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

  double get_energy_change(int r_d, int r_a, int index_d, int index_a,
                           double dc,
                           state_struct &state, 
                           interactions_struct &interactions){
    
    //std::cout << "Started get_energy_change\n";
    vec1i n_d = get_neighbours(r_d,state.Lx,state.Ly,state.Lz);
    vec1i n_a;

    bool nearest_neighbours = false;

    if(r_d==r_a){
      n_a = n_d;
    }
    else{
      n_a = get_neighbours(r_a,state.Lx,state.Ly,state.Lz);
      for(int i=0;i<n_d.size();i++){
        if(r_a==n_d[i]){
          nearest_neighbours = true;
        }
      }
    }

    double dE = 0;

    //std::cout << "Calculated neighbours\n";
    for(int j=0;j<n_d.size();j++){

      int r_d_j = n_d[j];
      int r_a_j = n_a[j];

      //std::cout << "r_d_j = " << r_d_j << "\n";
      //std::cout << "r_a_j = " << r_a_j << "\n";

      vec1d c_d_j = state.concentration[r_d_j];
      vec1d c_a_j = state.concentration[r_a_j];

      for(int b=0;b<interactions.coupling_matrix[j].size();b++){
        //std::cout << "b = " << b << "\n";
        dE-= interactions.coupling_matrix[j][index_d][b]*dc*c_d_j[b];
        dE+= interactions.coupling_matrix[j][index_a][b]*dc*c_a_j[b];
      }
    }
    //std::cout << "Starting nn correction\n";
    if(nearest_neighbours){
      int dr = get_bond_direction(r_d,r_a,state.Lx,state.Ly,state.Lz);
      dE-= interactions.coupling_matrix[dr][index_d][index_a]*dc*dc;
    }

    //std::cout << "Finished get_energy_change\n";
    return dE;
  }

  double get_entropy_change_shift(int r, int index, double dc,
                                  state_struct &state,
                                  interactions_struct &interactions,
                                  double eps){
    /*
     * Calculate the change in mixing entropy on a selected site as a result of
     * local density shift
     */

    double dS;

    // Define the concentration components before and after the shift
    double c_old = state.concentration[r][index];
    double c_new = c_old + dc;

    // Calculate the local density before and after the shift
    double rho_old = state.local_density[r];
    double rho_new = rho_old + dc;

    if(dc<0){
      // Calculate entropy change for the donor site
      dS = interactions.T_model*dc*log(c_old/(1-rho_old));

      if(c_new>=eps){
        dS+= interactions.T_model*c_new*log(c_new/c_old);
      }

      if(1-rho_new>=eps){
        dS+= interactions.T_model*(1-rho_new)*log((1-rho_new)/(1-rho_old));
      }
    }

    if(dc>0){
      // Calculate entropy change for the acceptor site
      dS = interactions.T_model*dc*log(c_new/(1-rho_new));

      if(c_old>=eps){
        dS+= interactions.T_model*c_old*log(c_new/c_old);
      }

      if(1-rho_old>=eps){
        dS+= interactions.T_model*(1-rho_old)*log((1-rho_new)/(1-rho_old));
      }
    }
    return dS;
  }

  double get_entropy_change_convert(int r, int index_d, int index_a, double dc,
                                  state_struct &state,
                                  interactions_struct &interactions,
                                  double eps){
    /*
     * Calculate the change in mixing entropy on a selected site as a result of
     * local re-distribution of fractional concentrations
     */

    double dS;

    // Define the donor concentration before and after the conversion
    double c_d_old = state.concentration[r][index_d];
    double c_d_new = c_d_old-dc;

    // Define the acceptor concentration before and after the conversion
    double c_a_old = state.concentration[r][index_a];
    double c_a_new = c_a_old+dc;

    if(c_a_old>=eps and c_d_new>=eps){
      // Generic concentration transfer
      dS = 0.5*interactions.T_model*dc*log(c_a_old*c_a_new/c_d_old/c_d_new); 
      dS+= 0.5*interactions.T_model*(c_a_old+c_a_new)*log(c_a_new/c_a_old);
      dS+= 0.5*interactions.T_model*(c_d_old+c_d_new)*log(c_d_new/c_d_old);
    }

    else if(c_a_old<eps and c_d_new>=eps){
      // Transfer part of concentration from donor to empty acceptor
      dS = interactions.T_model*dc*log(c_a_new/c_d_new); 
      dS+= interactions.T_model*c_d_old*log(c_d_new/c_d_old);
    }

    else if(c_a_old>=eps and c_d_new<eps){
      // Transfer everything from donor to non-empty acceptor
      dS = interactions.T_model*dc*log(c_a_old/c_d_old); 
      dS+= interactions.T_model*c_a_new*log(c_a_new/c_a_old);
    }

    else if(c_a_old<eps and c_d_new<eps){ 
      // Transfer everything from donor to empty acceptor
      dS = 0;
    }
    return dS;
  }

  void update_interactions(double dE, double dS, double dF, 
                           interactions_struct &interactions){
  
    interactions.energy += dE;
    interactions.entropy += dS;
    interactions.free_energy += dF;
  }
}
