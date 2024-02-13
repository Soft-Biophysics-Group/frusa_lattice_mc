#include "models.h" 

namespace simulation{
  /*
   * Definitions for the particle class
   */
  
  particles::particles(const model_data &model_data_1d,
                       const mc_data &mc_data_1d) :
    N(model_data_1d.N),
    Np(model_data_1d.Np),
    k11(model_data_1d.k11),
    k12(model_data_1d.k12),
    k21(model_data_1d.k21),
    rng(model_data_1d.rng)
    //T(mc_data_1d.T) 
  {
    /*
     * Initialize the system of particles and calculate the initial energy
     */


    /*Initialize pseudorandom number distributions*/
    uniform_dist  = real_dist(0,1);
    particle_dist = int_dist(0,Np-1);
    empty_dist    = int_dist(0,N-Np-1);
    binary_dist   = int_dist(0,1);

    initialize();
    coupling_matrix = {{{k11,k21},{k12,k11}},
                       {{k11,k12},{k21,k11}}};

    energy = get_energy();

}

  /*
   * Required public routines
   */

  void particles::print_state(){
    /*
     * Print the current state of the system
     */
    std::cout << "\n";
    for(int i=0;i<Np;i++){
      std::cout << positions[i] << " " << orientations[i] << "\n";
    }
    std::cout << "\n";
  }

  void particles::save_state(std::string file_name, std::string address){
    /*
     * Save the current state of the system to a file
     */
    std::ofstream state_f;
    state_f.open(address+file_name);
    if(!state_f){
      std::cerr << "Could not open "+address+file_name << std::endl;
      exit(1);
    }

    for(int i=0;i<Np;i++){
      state_f << positions[i] << " " << orientations[i] << "\n";
    }
    state_f.close();
  }
  
  void particles::print_energy(){
    /*
     * Print the current energy of the system
     */
    std::cout << energy << "\n";
  }

  void particles::update_state(double T){
    /*
     * Update the state of the system using Metropolis algorithm
     */
    for(int i=0;i<Np;i++){
      int particle_index = particle_dist(rng);
      int update_type = binary_dist(rng);
      
      if(update_type==0){
        update_position(particle_index,T);
      }
      else{
        update_orientation(particle_index,T);
      }
    } 
  }

  void particles::update_averages(){
  }

  void particles::save_averages(){
  }
  /*
   * Class-specific routines
   */

  void particles::initialize(){
    /*
     * Initialize a system of particles at random locations and with
     * random orientations 
     */
    
    /*Define a re-shuffled vector with all positions on a 1D lattice*/
    vec1i all_positions(N);

    for(int i=0;i<N;i++){
      all_positions[i] = i;
    }
    
    std::shuffle(std::begin(all_positions), std::end(all_positions), rng);

    /*Define and populate the state vector*/

    for(int i=0;i<Np;i++){
      
      positions.push_back(all_positions[i]);

      int orientation_type = binary_dist(rng);

      if(orientation_type==0){
        orientations.push_back(1);
      }
      else if(orientation_type==1){
        orientations.push_back(2);
      }
    }
  }

  vec2i particles::get_neighbours(int r){
    
    /*
     * Extract the indices (location in the state vector) of the neighbours
     * of a specified particle at position r
     */

    vec2i neighbours;

    int rm = (N+((r-1)%N))%N;
    int rp = (r+1)%N;

    for(int i=0;i<Np;i++){
      if(positions[i]==rm){
        neighbours.push_back({0,i});
      }
      else if(positions[i]==rp){
        neighbours.push_back({1,i});
      }
    }

    return neighbours;
  }

  double particles::get_energy(){
    /*
     * Calculate total energy of the system
     */

    double en = 0;

    for(int i=0;i<Np;i++){
      
      vec2i n = get_neighbours(positions[i]);

      for(int j=0;j<n.size();j++){
        int k = n[j][0];
        en+= coupling_matrix[k][ orientations[i]-1][ orientations[n[j][1]]-1];
      }
    }
    return en/2;
  }

  void particles::update_position(int particle_index, double T){
    /*
     * Attempt to move the selected particle (particle_index) to a new position
     *  with Metropolis acceptance rate
     */
    vec1i empty_sites;

    for(int i=0;i<N;i++){
      if(std::find(positions.begin(), positions.end(), i) == positions.end()){
        /* Position vector doesn't contain i */
        empty_sites.push_back(i);
      }
    }

    int empty_index = empty_dist(rng);
    
    int position_old = positions[particle_index];
    int position_new = empty_sites[empty_index];

    /*Find the neighbours of the old and new positions*/
    vec2i n_old = get_neighbours(position_old);
    vec2i n_new = get_neighbours(position_new);

    /*Check if the old and new positions are first neighbours*/
    bool direct_neighbours[2];

    int rp1 = (position_old+1)%N;
    int rm1 = (N+((position_old-1)%N))%N;

    if(position_new==rp1){
     direct_neighbours[0] = false;
     direct_neighbours[1] = true;
    }
    else if(position_new==rm1){
      direct_neighbours[0] = true;
      direct_neighbours[1] = false;
    }
    else{
      direct_neighbours[0] = false;
      direct_neighbours[1] = false;
    }

    /*Calculate the energy difference resulting from the position change*/
    double dE = 0;
    
    for(int i=0; i<n_new.size();i++){
      if(not direct_neighbours[i]){
        /*Exclude the possibility of interaction with itself*/
        int k = n_new[i][0];
        dE+= coupling_matrix[k][ orientations[particle_index]-1 ]\
                               [ orientations[n_new[i][1]]-1 ];  
      }
    }
    for(int i=0; i<n_old.size();i++){
      int k = n_old[i][0];
      dE-= coupling_matrix[k][ orientations[particle_index]-1 ]\
                               [ orientations[n_old[i][1]]-1 ];
    }

    /*Accept the new position using Metropolis rule*/
    if(dE<=0){
      positions[particle_index] = position_new;
      energy+=dE;
    }
    else if(exp(-dE/T)>uniform_dist(rng)){
      positions[particle_index] = position_new;
      energy+=dE;
    }
  }

  void particles::update_orientation(int particle_index, double T){
    /* 
     * Attempt to change the orientation of the selected particle 
     * (particle_index) with Metropolis acceptance rate
     */
    
    /*Determine the new proposed orientation*/
    int orientation_old = orientations[particle_index];
    int orientation_new;
    if(orientation_old==1){
      orientation_new=2;
    }
    else{
      orientation_new=1;
    }

    /*Calculate the energy cost of changing the orientation*/
    double dE = 0;

    vec2i n = get_neighbours(positions[particle_index]);

    for(int i=0; i<n.size();i++){
      int k = n[i][0];
      dE+= coupling_matrix[k][ orientation_new-1]\
                             [ orientations[n[i][1]]-1];
      dE-= coupling_matrix[k][ orientation_old-1]\
                             [ orientations[n[i][1]]-1];
    }

    /*Accept the new position using Metropolis rule*/
    if(dE<=0){
      orientations[particle_index] = orientation_new;
      energy+=dE;
    }
    else if(exp(-dE/T)>uniform_dist(rng)){
      orientations[particle_index] = orientation_new;
      energy+=dE;
    }


  }

  void particles::update_psi(vec2d &psi){
    /*
     * Update the density vector
     */

    for(int i=0;i<Np;i++){
      psi[positions[i]][orientations[i]-1] += 1;
    }
  }
}

