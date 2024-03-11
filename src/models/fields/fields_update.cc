#include "fields_update.h"

namespace fields_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
  void update_system(state_struct &state, 
                     interactions_struct &interactions,
                     model_parameters_struct &parameters,
                     double T){
    
    int_dist binary_dist(0,1);

    for(int r=0;r<state.N;r++){
      int update_type = binary_dist(parameters.rng);

      if(update_type==0){
        shift_local_density(state,interactions,parameters,T);
      }
      else{
        convert_concentrations(state,interactions,parameters,T);
      }
    } 
  }
  /* 
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
  void shift_local_density(state_struct &state, 
                           interactions_struct &interactions,
                           model_parameters_struct &parameters,
                           double T, double eps){
    
    real_dist uniform_dist(0,1);

    int_dist donor_dist(0,state.donor_list.size()-1);
    int_dist acceptor_dist(0,state.acceptor_list.size()-1);

    int d_list_ind = donor_dist(parameters.rng);
    int a_list_ind = acceptor_dist(parameters.rng);

    int r_d = state.donor_list[d_list_ind];
    int r_a = state.acceptor_list[a_list_ind];
   
    while(r_d==r_a){
      r_a = state.acceptor_list[acceptor_dist(parameters.rng)];
    }
    
    int index_d = select_element(r_d,0,state,parameters);
    int index_a = select_element(r_a,1,state,parameters);

    double c_d = state.concentration[r_d][index_d];
    double c_a = state.concentration[r_a][index_a];
    
    double bound_d = c_d;
    double bound_a = std::min(1-c_a,1-state.local_density[r_a]);

    double bound_total = std::min(bound_d,bound_a);

    double dc;

    if(bound_total<eps){
      // No need to shift tiny amounts of density (below threshold eps), 
      // simply transfer all density from the site
      dc = bound_total;
    }
    else{
      dc = uniform_dist(parameters.rng)*bound_total;
    }

    double dE = get_energy_change(r_d,r_a,index_d,index_a,dc,
                                  state,interactions);

    double dS = 0;
    dS += get_entropy_change_shift(r_d,index_d,-dc,state,interactions);
    dS += get_entropy_change_shift(r_a,index_a, dc,state,interactions);

    double dF = dE+dS;

    // Accept the new position using Metropolis rule
    if(dF<=0){

      update_state(r_d,index_d,d_list_ind,-dc,state);
      update_state(r_a,index_a,a_list_ind, dc,state);

      update_interactions(dE,dS,dF,interactions);
    }
    else if(exp(-dF/T)>uniform_dist(parameters.rng)){

      update_state(r_d,index_d,d_list_ind,-dc,state);
      update_state(r_a,index_a,a_list_ind, dc,state);

      update_interactions(dE,dS,dF,interactions);
    }
  }
  
  void convert_concentrations(state_struct &state, 
                              interactions_struct &interactions,
                              model_parameters_struct &parameters,
                              double T, double eps){

    real_dist uniform_dist(0,1);

    int_dist donor_dist(0,state.donor_list.size()-1);

    int list_ind = donor_dist(parameters.rng);

    int r = state.donor_list[list_ind];
   
    int index_d = select_element(r,0,state,parameters);
    int index_a = select_element(r,1,state,parameters);

    while(index_d==index_a){
      index_a = select_element(r,1,state,parameters);
    }

    double c_d = state.concentration[r][index_d];
    double c_a = state.concentration[r][index_a];
    
    double bound_d = c_d;
    
    double dc;

    if(bound_d<eps){
      // No need to shift tiny amounts of density (below threshold eps), 
      // simply convert all concentration to the other type
      dc = bound_d;
    }
    else{
      dc = uniform_dist(parameters.rng)*bound_d;
    }
    
    double dE = get_energy_change(r,r,index_d,index_a,dc,
                                  state,interactions);

    double dS = 0;
    dS += get_entropy_change_convert(r,index_d,index_a,dc,state,interactions);

    double dF = dE+dS;

    // Accept the new position using Metropolis rule
    if(dF<=0){

      update_state(r,index_d,list_ind,-dc,state);
      update_state(r,index_a,list_ind, dc,state);

      update_interactions(dE,dS,dF,interactions);
    }
    else if(exp(-dF/T)>uniform_dist(parameters.rng)){

      update_state(r,index_d,list_ind,-dc,state);
      update_state(r,index_a,list_ind, dc,state);

      update_interactions(dE,dS,dF,interactions);
    }
  }

  int select_element(int r, double bound,
                     state_struct &state, 
                     model_parameters_struct &parameters){
    
    vec1d c = state.concentration[r];
    vec1i c_ind;

    for(int s=0;s<state.ns;s++){
      if(c[s]!=bound){
        c_ind.push_back(s);
      }
    }
   
    int_dist c_dist(0,c.size()-1);
    
    return c_ind[c_dist(parameters.rng)];
  }
}
