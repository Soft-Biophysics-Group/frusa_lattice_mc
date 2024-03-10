namespace fields_space{
  
  /*
   * Definitions required for the public routines of the model class
   */
 void update_state(); 
  /* 
   * End of the required definitions for the model class
   */

  /*
   * Library-specific definitions
   */
  void shift_local_density(state_struct &state, 
                           interactions_struct &interactions,
                           model_parameters_struct &parameters,
                           double T){
    
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
    double bound_a = std::min(1-c_a;1-state.local_density[r_a]);

    bound_total = std::min(bound_d,bound_a);

    if(bound_total<eps){
      // No need to shift tiny amounts of density (below threshold eps), 
      // simply transfer all density from the site
      drho = bound_total;
    }
    else{
      drho = uniform_dist(parameters.rng)*bound_total;
    }

    vec1i n_d = get_neighbours(r_d);
    vec1i n_a = get_neighbours(r_a);

    // Check if the donor and the acceptor are nearest-neighbours
    bool nearest_neighbours;

    for(int i=0;i<n_d.size();i++){
      if(r_a==n_d[i]){
        nearest_neighbours = true;
      }
    }

    double dE = 0;
    double dS = 0;

    for(int j=0;j<n_d.size();j++){
      
      int r_d_j = n_d[j];
      int r_a_j = n_a[j];
      
      vec1d c_d_j = state.concentration[r_d_j];
      vec1d c_a_j = state.concentration[r_a_j];

      for(int b=0;b<state.ns;b++){
        dE+= interactions.coupling_matrix[j][index_a][b]*drho*c_a_j[k];
        dE-= interactions.coupling_matrix[j][index_d][b]*drho*c_d_j[k];
      }
    }
    if(nearest_neighbours){
      dr = get_bond_direction(r_d,r_a);
      dE-= interactions.coupling_matrix[dr][index_d][index_a]*drho*drho;
    }

    dS += get_entropy_shift(c_d,index_d,-drho);
    dS += get_entropy_shift(c_a,index_a, drho);

    double dF = dE+dS;

    /*Accept the new position using Metropolis rule*/
    if(dF<=0){

      state.concentration[r_d][index_d] -= drho;
      state.concentration[r_a][index_a] += drho;
      
      state.local_density[r_d] -= drho;
      state.local_density[r_a] += drho;

      if(state.local_density[r_d]==0){
        state.donor_list.erase(d_list_ind);
      }
      if(state.local_density[r_a]==1){
        state.acceptor_list.erase(a_list_ind);
      }
      
      interactions.energy      += dE;
      interactions.entropy     += dS;
      interactions.free_energy += dF;
    
    }
    else if(exp(-dF/T)>uniform_dist(parameters.rng)){

      state.concentration[r_d][index_d] -= drho;
      state.concentration[r_a][index_a] += drho;
      
      state.local_density[r_d] -= drho;
      state.local_density[r_a] += drho;

      if(state.local_density[r_d]==0){
        state.donor_list.erase(d_list_ind);
      }
      if(state.local_density[r_a]==1){
        state.acceptor_list.erase(a_list_ind);
      }
      
      interactions.energy      += dE;
      interactions.entropy     += dS;
      interactions.free_energy += dF;
    
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
    
    return c_ind[c_dist(params.rng)];
  }
}
